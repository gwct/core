#!/usr/bin/env python3

"""
Extract all iMessage/SMS attachments from an iOS backup created with:

    idevicebackup2 backup /mnt/nell-hdd/Archive/iphone-backups/<YYYYMM>

This script assumes:

  - You choose which month to use (BACKUP_TIMESTAMP, e.g. "202512").
  - Under that month directory, idevicebackup2 created one or more
    "UDID" subdirectories (like 00008030-...).
  - Each UDID subdir contains a Manifest.db.

The script:

  - Finds the newest UDID subdir under the chosen month (by Manifest.db mtime).
  - Locates the Messages DB (sms.db) via Manifest.db.
  - Reads attachment entries from sms.db.
  - Copies all attachment files into a single output dir:

        /mnt/nell-hdd/Archive/iphone-backups/<YYYYMM>/iMessage_attachments
"""

import os
import shutil
import sqlite3
from pathlib import Path
from typing import Optional, Tuple

# ---------------------------------------------------------------------------
# Configuration (you change only BACKUP_TIMESTAMP)
# ---------------------------------------------------------------------------

# Base directory under which you store month-stamped iPhone backups.
BACKUP_BASE = Path("/mnt/nell-hdd/Archive/iphone-backups")

# The month you want to use (YYYYMM string).
# You can change this per run, e.g. "202512", "202601", etc.
BACKUP_TIMESTAMP = "202512"

# Name of the output directory inside the chosen month directory
# where extracted attachments will be written.
OUTPUT_DIR_NAME = "iMessage-attachments"


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def find_latest_udid_backup_dir(month_dir: Path) -> Optional[Path]:
    """
    Given a month directory (e.g. /.../iphone-backups/202512),
    find the UDID backup subdirectory that has the newest Manifest.db.

    Example structure:
        /.../iphone-backups/202512/
            00008030-001A69E001C3402E/
                Manifest.db
                ...
            another-udid-.../
                Manifest.db
                ...

    Returns:
        backup_dir (Path) or None if no Manifest.db found under month_dir.
    """
    manifest_paths = []

    # Only look one level down (UDID dirs) rather than walking entire tree:
    # if you might have nested stuff, you can switch to os.walk().
    for child in month_dir.iterdir():
        if not child.is_dir():
            continue
        m_path = child / "Manifest.db"
        if not m_path.exists():
            continue
        try:
            mtime = m_path.stat().st_mtime
        except OSError:
            continue
        manifest_paths.append((mtime, m_path))

    if not manifest_paths:
        return None

    # Pick the Manifest.db with newest modification time
    manifest_paths.sort(key=lambda x: x[0], reverse=True)
    _, newest_manifest = manifest_paths[0]

    backup_dir = newest_manifest.parent  # the UDID folder
    return backup_dir


def get_file_id_by_relpath(cur, logical_path: str) -> Optional[str]:
    """
    Look up a fileID by relativePath only (used for sms.db, which we already know).

    For example:
      logical_path = 'Library/SMS/sms.db'
    """
    cur.execute(
        "SELECT fileID FROM Files WHERE relativePath = ?",
        (logical_path,)
    )
    row = cur.fetchone()
    return row[0] if row else None


def get_file_id_by_domain_and_relpath(cur, domain: str, logical_path: str) -> Optional[str]:
    """
    Look up a fileID by (domain, relativePath), e.g.:

      domain = 'MediaDomain'
      logical_path = 'Library/SMS/Attachments/xx/yy/...'

    Returns the fileID string or None if not found.
    """
    cur.execute(
        "SELECT fileID FROM Files WHERE domain = ? AND relativePath = ?",
        (domain, logical_path)
    )
    row = cur.fetchone()
    return row[0] if row else None


# ---------------------------------------------------------------------------
# Main extraction logic
# ---------------------------------------------------------------------------

def main():
    # 1. Resolve the month directory from BACKUP_TIMESTAMP
    month_dir = BACKUP_BASE / BACKUP_TIMESTAMP

    if not month_dir.exists() or not month_dir.is_dir():
        print(f"ERROR: Month directory does not exist or is not a directory: {month_dir}")
        return

    # 2. Find the latest UDID backup directory under that month
    backup_dir = find_latest_udid_backup_dir(month_dir)
    if backup_dir is None:
        print(f"ERROR: No Manifest.db found under {month_dir}")
        return

    manifest_db = backup_dir / "Manifest.db"
    output_dir = month_dir / OUTPUT_DIR_NAME

    if not manifest_db.exists():
        print("ERROR: Manifest.db not found in backup dir:", manifest_db)
        return

    print("Using month dir :", month_dir)
    print("Using backup dir:", backup_dir)
    print("Manifest DB     :", manifest_db)
    print("Output dir      :", output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 3. Open Manifest.db
    conn = sqlite3.connect(str(manifest_db))
    cur = conn.cursor()

    # 4. Find the Messages database (sms.db) in the backup
    messages_logical_path = "Library/SMS/sms.db"
    sms_db_id = get_file_id_by_relpath(cur, messages_logical_path)
    if not sms_db_id:
        print("ERROR: Could not find sms.db in Manifest.db – is this a valid backup?")
        return

    # In the backup, the file is stored under:
    #   <backup_dir>/<first 2 chars of fileID>/<fileID>
    sms_db_path = backup_dir / sms_db_id[:2] / sms_db_id
    if not sms_db_path.exists():
        print("ERROR: sms.db file not found at expected path:", sms_db_path)
        return

    print("Found Messages DB at:", sms_db_path)

    # 5. Open the Messages database
    sms_conn = sqlite3.connect(str(sms_db_path))
    sms_cur = sms_conn.cursor()

    # 6. Get all attachments that have a filename.
    #    'attachment' table has columns including:
    #        ROWID, filename, mime_type
    #    filename is like '~/Library/SMS/Attachments/...'
    sms_cur.execute("""
        SELECT ROWID, filename, mime_type
        FROM attachment
        WHERE filename IS NOT NULL
    """)
    rows = sms_cur.fetchall()
    print(f"Found {len(rows)} attachment records in sms.db.")

    copied = 0
    missing = 0

    # 7. For each attachment, locate the actual file in the backup and copy it
    for rowid, filename, mime in rows:
        if not filename:
            continue

        # filename is usually like:
        #   '~/Library/SMS/Attachments/77/07/BE83.../image000000.jpg'
        #
        # In Manifest.db, relativePath is stored as:
        #   'Library/SMS/Attachments/77/07/BE83.../image000000.jpg'
        #
        # So we strip '~/'
        path = filename

        # Strip leading "~/"
        if path.startswith("~/"):
            path = path[2:]
        # Also strip any leading "/" just in case
        path = path.lstrip("/")

        # This should now match the Files.relativePath
        logical_path = path  # e.g. 'Library/SMS/Attachments/...'

        # Attachments live under domain = 'MediaDomain'
        file_id = get_file_id_by_domain_and_relpath(cur, "MediaDomain", logical_path)
        if not file_id:
            # Not all attachments will necessarily be present in the backup
            missing += 1
            continue

        # Backup storage path: <backup_dir>/<first 2 chars of fileID>/<fileID>
        src = backup_dir / file_id[:2] / file_id
        if not src.exists():
            missing += 1
            continue

        # Build an output filename:
        #  - Use ROWID to avoid collisions
        #  - Preserve original extension (if any)
        _, ext = os.path.splitext(filename)
        out_name = f"attachment_{rowid}{ext}"
        dst = output_dir / out_name

        try:
            shutil.copy2(src, dst)
            copied += 1
        except Exception as e:
            print(f"Error copying {src} -> {dst}: {e}")

    print(f"Copied {copied} attachments to {output_dir}")
    print(f"Skipped/missing: {missing}")


if __name__ == "__main__":
    main()