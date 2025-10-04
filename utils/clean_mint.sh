#!/bin/bash
# =======================================
# Script to clean up old kernels and other unnecessary files on a Debian-based system.
# It keeps the current and previous kernels, and prompts the user before deleting anything.
# Usage: Run the script in a terminal. It may require sudo privileges for some operations.
# Note: Uncomment the thumbnail cache section if you want to enable that cleanup.
# Date: 2025-10
# =======================================

set -e

echo "Before cleanup:"
root_info=$(df -h / | awk 'NR==2 {print $4 " free / " $2 " total (" $5 " used)"}')
echo "Root filesystem: $root_info"

# =======================================
# ========== Kernel cleanup =============

echo "Detecting installed kernels..."
# Get list of installed kernels, sorted by version (oldest first)
installed_kernels=($(dpkg --list | awk '/^ii/ && $2 ~ /^linux-image-[0-9]/ { print $2 }' | sort -V))

# Get running kernel version (as package name)
current_kernel="linux-image-$(uname -r)"

# Remove current kernel from the list
filtered_kernels=()
for k in "${installed_kernels[@]}"; do
    if [[ "$k" != "$current_kernel" ]]; then
        filtered_kernels+=("$k")
    fi
done
# Find previous kernel: the highest version kernel not running
if [ "${#filtered_kernels[@]}" -gt 0 ]; then
    previous_kernel="${filtered_kernels[-1]}"
else
    previous_kernel=""
fi

# List for removal: all except current and previous
removal_kernels=()
for k in "${installed_kernels[@]}"; do
    if [[ "$k" != "$current_kernel" && "$k" != "$previous_kernel" ]]; then
        removal_kernels+=("$k")
    fi
done

if [ "${#removal_kernels[@]}" -eq 0 ]; then
    echo "No old kernels to remove."
else
    echo "Current kernel: $current_kernel"
    echo "Previous kernel: $previous_kernel"
    echo "The following old kernel(s) will be removed:"
    for k in "${removal_kernels[@]}"; do
        echo "  $k"
    done

    # Estimate space to be recovered
    kernel_files=()
    for pkg in "${removal_kernels[@]}"; do
        files=$(dpkg -L "$pkg" 2>/dev/null | grep -E 'vmlinuz|initrd.img|System.map|abi|modules$')
        kernel_files+=($files)
    done

    bytes=0
    for f in "${kernel_files[@]}"; do
        if [ -e "$f" ]; then
            b=$(stat -c%s "$f" 2>/dev/null || echo 0)
            bytes=$((bytes+b))
        fi
    done
    mb=$((bytes/1024/1024))
    echo "Estimated space to be freed (kernel files): ${mb} MB"

    read -p "Proceed with removal of these kernels? [y/N] " confirm
    if [[ "$confirm" =~ ^[Yy]$ ]]; then
        sudo apt remove --purge "${removal_kernels[@]}"
        sudo apt autoremove --purge
        echo "Old kernels removed."
    else
        echo "Skipped kernel removal."
    fi
fi

# ========== Other cleaning prompts ==========

# 1. Clean APT cache
apt_cache_size=$(sudo du -sh /var/cache/apt/archives | cut -f1)
echo
echo "APT package cache location: /var/cache/apt/archives"
echo "Current size: $apt_cache_size"
read -p "Clean APT cache to free up space? [y/N] " clean_apt
if [[ "$clean_apt" =~ ^[Yy]$ ]]; then
    sudo apt clean
    echo "APT cache cleaned."
fi

# 2. Clean thumbnail cache
# thumb_cache="$HOME/.cache/thumbnails"
# if [ -d "$thumb_cache" ]; then
#     thumb_size=$(du -sh "$thumb_cache" | cut -f1)
#     echo
#     echo "Thumbnail cache: $thumb_cache"
#     echo "Current size: $thumb_size"
#     read -p "Clean thumbnail cache? [y/N] " clean_thumbs
#     if [[ "$clean_thumbs" =~ ^[Yy]$ ]]; then
#         rm -rf "$thumb_cache"/*
#         echo "Thumbnail cache cleaned."
#     fi
# fi

# 3. Clean journal logs
if [ -d /var/log/journal ]; then
    journal_size=$(sudo du -sh /var/log/journal | cut -f1)
    echo
    echo "Systemd journal logs location: /var/log/journal"
    echo "Current size: $journal_size"
    read -p "Clean journal logs (keep 100M)? [y/N] " clean_journal
    if [[ "$clean_journal" =~ ^[Yy]$ ]]; then
        sudo journalctl --vacuum-size=100M
        echo "Journal logs cleaned (retaining 100M)."
    fi
fi

# =======================================

echo
echo "After cleanup:"
root_info=$(df -h / | awk 'NR==2 {print $4 " free / " $2 " total (" $5 " used)"}')
echo "Root filesystem: $root_info"

echo
echo "Cleanup complete."