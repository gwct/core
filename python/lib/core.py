#############################################################################
# CORE functions
# Gregg Thomas
# August 2013-present
#############################################################################

import string, sys, os, re, subprocess, datetime, gzip, random
from collections import defaultdict
import timeit

#############################################################################

def loadingBar(counter, length, done, bars, firstbar=False, disperc=False):
# This function serves as a text loading bar for long scripts with counters. The following
# lines must be added within the script to initialize and terminate the script:
# Initilization:
# numlines = core.getFileLen(alnfilename);
# numbars = 0;
# donepercent = [];
# i = 0;
# Termination:
#    pstring = "100.0% complete.";
#    sys.stderr.write('\b' * len(pstring) + pstring);
#    print "\nDone!";
#
# If length is lines in a file use the core.getFileLen function to get the number of lines in the file

    # try:
    #     if sys.version[0] == '2':
    #         pchr = u'\u2591'.encode('utf-8');
    #         lchr = u'\u2588'.encode('utf-8');
    #     elif sys.version[0] == '3':
    #         pchr = u'\u2591';
    #         lchr = u'\u2588';        
    # except:
    #     pchr, lchr = "*","*";

    try:
        #pchr, lchr=u'\u2591',u'\u2588';
        pchr, lchr=u'\u2588',u'\u2591';
    except:
        pchr, lchr = "=","|";

    percent = float(counter) / float(length) * 100.0;
    percentdone = int(percent);

    p = str(percent)
    pstring = " " + p[:5] + "% complete.";

    if percentdone % 2 == 0 and done != None and percentdone not in done:
        loading = "|";
        j = 0;
        while j < bars:
            loading += pchr;
            j += 1;
        if j <= 49:
            loading += lchr;
        else:
            loading += pchr;
        j += 1;
        if j == 50:
            loading = loading[:-1] + pchr;

        while j < 50:
            loading += "-";
            j += 1;
        loading += "|";

        if disperc:
            loading += "                 ";
        if firstbar:
            sys.stderr.write(loading);
            firstbar = False
        else:
            sys.stderr.write('\b' * len(loading) + loading);

        done.append(percentdone);
        bars = bars + 1;
    if disperc:
        sys.stderr.write('\b' * len(pstring) + pstring);
    sys.stderr.flush();
    
    return bars, done, firstbar;

#############################################################################

def loadingRotator(counter, rotate, divisor):
# Provides a loading rotator for loops. The following line must be used to initialize the function
# before the loop in the main code:
# rotator = 0;

    rotation = ['|', '/', '-', '\\'];

    if counter % divisor == 0:
        sys.stderr.write('\b' + rotation[rotate]);
        rotate = rotate + 1;
        if rotate >= len(rotation):
            rotate = 0;

    return rotate;

#############################################################################

def mean(data):
# Calculates and returns the mean of a list of numbers.
    return sum(data) / len(data);

#############################################################################

def variance(data):
# Calculates and returns the variance of a list of numbers.
    mean = mean(data);
    var = 0.0;
    for d in data:
        var = var + (float(d) - mean)**2;
    var = var / ((float(len(data)) - 1.0));
    return var;

#############################################################################

def getFileLen(i_name):
# Calls 'wc -l' to get the number of lines in the file.
    p = subprocess.Popen(['wc', '-l', i_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE);
    result, err = p.communicate();
    if p.returncode != 0:
        raise IOError(err);
    return int(result.strip().split()[0]);

#############################################################################

def getFileLenRead(i_name):
# Gets the number of lines in a file.
    num_lines = 0;
    for line in getFileReader(i_name)(i_name): num_lines += 1;
    return float(num_lines);

#############################################################################

def getDateTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m.%d.%Y | %H:%M:%S");

#############################################################################

def getDate():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m.%d.%Y");

#############################################################################

def getTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%I:%M:%S");

#############################################################################

def getLogTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m.%d.%Y-%I.%M.%S");

#############################################################################

def printWrite(o_name, o_line, file_flag=True):
# Function to print a string AND write it to the file.
    print(o_line);
    if file_flag == False:
        with open(o_name, "a") as f:
            f.write(o_line + "\n");

#############################################################################

def PW(o_line, o_name, file_flag=True):
# Function to print a string AND write it to the file. Redundant with printWrite
# but with a shorter name
    print(o_line);
    if file_flag == True:
        with open(o_name, "a", encoding="utf-8") as f:
            f.write(o_line + "\n");

#############################################################################

def PWS(o_line, o_stream=False, std_stream=True):
# Function to print a string AND write it to the file. More functional than above.
    if std_stream:
        print(o_line);
    if o_stream:
        o_stream.write(o_line + "\n");

#############################################################################

def logCheck(lopt, lfilename, outline):
# Function checks whether or not to write to a logfile, print something, or both.
    if lopt == 1:
        printWrite(lfilename, outline);
    else:
        print(outline);

#############################################################################

def errorOut(errnum, errmsg, ropt=0):
# Formatting for error messages.
    fullmsg = "| ** Error " + str(errnum) + ": " + errmsg + " |";
    border = " " + "-" * (len(fullmsg)-2);
    if ropt:
        return "\n" + border + "\n" + fullmsg + "\n" + border + "\n";
    else:
        print("\n" + border + "\n" + fullmsg + "\n" + border + "\n");

#############################################################################

def getOutdir(indir, prefix, suffix, stime):
# Retrieves full input directory name and proper output directory name for other scripts.
    if not os.path.isdir(indir):
        errorOut(0, "-i must be a valid directory path");
        sys.exit();
    indir = os.path.abspath(indir);
    filelist = os.listdir(indir);
    used = [0];
    for each in filelist:
        if each.find("-" + prefix) != -1:
            used.append(int(each[:each.index("-")]));
    outdir = os.path.join(indir, str(max(used)+1) + "-" + prefix + "-" + stime + suffix);

    return indir, outdir;

#############################################################################

def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a string to make it a given length
    spaces = sep * (totlen - len(string));
    return string + spaces;

#############################################################################

def filePrep(filename, header=""):
# Writes over a file, header optional (if no header just pass "")
    f = open(filename, "w");
    f.write(header);
    f.close();

#############################################################################

def listCheck(lst):
# Checks if all elements in a list are the same or not.
    if lst.count(lst[0]) == len(lst):
        return True;
    else:
        return False;

#############################################################################

def defaultOutFile(input_name, file_flag, suffix="", output_init=False):
    i = 2;
    if suffix != "" and suffix[0] != "-":
        suffix = "-" + suffix;
    if not output_init:
        output, ext = list(os.path.splitext(input_name));
    # If the user did not specify an output file name, take the base of the input file name.
    else:
        output, ext = list(os.path.splitext(output_init));
    # Otherwise, use the user specified option.
    if not file_flag:
        if output[-1] in ["\\", "/"]:
            output = output[:-1];
        output = output + suffix + "-1" + ".txt";
    else:
        output = output + suffix + "-1" + ext;

    while os.path.exists(output) or os.path.exists(os.path.splitext(output)[0]):
        output = os.path.splitext(output);
        output = output[0][:output[0].rindex("-")+1] + str(i) + output[1];
        i += 1;
    # If the chosen output file exists, this will continually add 1 to a counter label at the end
    # of the file until a new file is chosen that does not exist.

    return output, (i-1);

#############################################################################

def defaultOutDir(input_name, file_flag, suffix="", output_init=False):
    i = 2;
    if suffix != "" and suffix[0] != "-":
        suffix = "-" + suffix;
    if not output_init:
        output = os.path.splitext(input_name)[0].rstrip("/").rstrip("\\") + suffix;
    # If the user did not specify an output name, a directory will be made based on the input directory name.
    else:
        output = output_init;
    # Otherwise, use the user specified option.
    output += "-1";

    while os.path.exists(output):
        output = output[:output.rindex("-")+1] + str(i);
        i += 1;
    # If the chosen output directory exists, this will continually add 1 to a counter label at the end
    # of the directory until a new directory is chosen that does not exist.

    return output, (i-1);

#############################################################################

def getFileReader(i_name):
# Check if a file is gzipped, and if so set gzip as the file reader. Otherwise, read as a normal text file.
    try:
        gzip_check = gzip.open(i_name).read(1);
        reader = gzip.open;
    except:
        reader = open;
    return reader;

#############################################################################

def detectCompression(filename):
# Detect compression of a file by examining the first lines in the file

    compression_type = "none";

    magic_dict = {
            b"\x1f\x8b\x08": "gz",
            # b"\x1f\x8b\x08\x08": "gz",
            b"\x42\x5a\x68": "bz2",
            b"\x50\x4b\x03\x04": "zip"
        }
    # An encoded set of possible "magic strings" that start different types of compressed files
    # From: https://www.garykessler.net/library/file_sigs.html
    # \x is the escape code for hex values
    # b converts strings to bytes

    max_len = max(len(x) for x in magic_dict)
    # The number of characters to read from the beginning of the file should be the length of
    # the longest magic string

    file_start = open(filename, "rb").read(max_len);
    # Read the beginning of the file up to the length of the longest magic string

    for magic_string in magic_dict:
        if file_start.startswith(magic_string):
            compression_type = magic_dict[magic_string];
    # Check each magic string against the start of the file

    return compression_type;

#############################################################################
# The two line reader functions based on whether the input file is gzipped or not
# Better handled as lambda functions whereve they're needed.
def readLine(line):
    return line.strip().split("\t");

def readGzipLine(line):
    return line.decode().strip().split("\t");

#############################################################################

def dsum(*dicts):
# Given a list of dictionaries, this function merges them into a single dictionary, summing values of common keys.
    ret = defaultdict(int);
    for d in dicts:
        for k, v in d.items():
            ret[k] += v;
    return dict(ret);

#############################################################################

def chunks(l, n):
# Splits a list l into even chunks of size n.
    n = max(1, n)
    return (l[i:i+n] for i in range(0, len(l), n))

#############################################################################

def runTime(msg=False, writeout=False, printout=True):
# Some info to print at the beginning of a script
    if msg:
        if not msg.startswith("#"):
            msg = "# " + msg;
        PWS(msg, writeout, printout);

    PWS("# PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])), writeout, printout)
    PWS("# Script call:    " + " ".join(sys.argv), writeout, printout)
    PWS("# Runtime:        " + datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S"), writeout, printout);
    PWS("# ----------------", writeout, printout);

#############################################################################

def getRandStr(strlen=6):
# This function generates a random string to add onto the end of tmp files to avoid possible overwrites.
    return ''.join(random.choice(string.ascii_letters) for m in range(strlen));

#############################################################################

def report_step(step, pids, prog_start_time, step_start_time, step_status, logfilename=False, v=1, start=False, psutil_flag=False, full_update=False):
# Uses psutil to gather memory and time info between steps and print them to the screen.

    dashes = 150
    if psutil_flag:
        import psutil;
        dashes = 175;
    # Determine the number of dashes to frame the update table depending on the presence of psutil

    cur_time = timeit.default_timer();
    # The time at the start of the status update

    col_widths = [ 14, 10, 40, 40, 20, 16 ];
    if psutil_flag:
        col_widths += [25, 20];
    # The column widths

    if start:
        headers = [ "# Date", "Time", "Current step", "Status", "Elapsed time (s)", "Step time (s)" ];
        if psutil_flag:
            headers += ["Current mem usage (MB)", "Virtual mem usage (MB)"]
        # A list of the headers

        headers = "".join([ spacedOut(str(headers[i]), col_widths[i]) for i in range(len(headers)) ]);
        # Converting the list to a string based on the column widths

        if logfilename:
            printWrite(logfilename, v, "# " + "-" * dashes);
            printWrite(logfilename, v, headers);
            printWrite(logfilename, v, "# " + "-" * dashes);
        # Print the dashes and the headers
    # The first call is just to print the headers

    ##########

    else:
        prog_elapsed = str(round(cur_time - prog_start_time, 5));
        # Get the total amount of time that the program has been running

        if not step_start_time:
        # If no step start time is given, then this is the first entry for this status
        # update, that will display "In progress..." or similar.

            out_line = [ "# " + getDate(), getTime(), step, step_status ];
            # The output for the initial status entry includes the date, time, step label, and progress message

            term_col_widths = col_widths[:4];
            # Only get the first 4 column widths for the initial status entry

            out_line = [ spacedOut(str(out_line[i]), term_col_widths[i]) for i in range(len(out_line)) ];

            if full_update:
                out_line += "\n";
            # For some status updates, intermediate info will be printed, in which case we add a newline here

            sys.stdout.write("".join(out_line));
            sys.stdout.flush();
            # Convert the output list to a string, write, and flush stdout

        # The initial status entry to display "In progress..."

        #####

        else:
            step_elapsed = str(round(cur_time - step_start_time, 5));
            # Get the full step time here

            out_line = [ step_status, prog_elapsed, step_elapsed ];
            # Gather info for the full output line to print to screen

            if psutil_flag:
                mem = round(sum([p.memory_info()[0] for p in pids]) / float(2 ** 20), 5);
                vmem = round(sum([p.memory_info()[1] for p in pids]) / float(2 ** 20), 5);
                out_line += [str(mem), str(vmem)];
            # If psutil is present, get current memory info

            term_col_widths = col_widths[3:];
            # Get the column widths for the print to screen output

            file_line = [ "# " + getDate(), getTime(), step ] + out_line;
            file_col_widths = col_widths[:3] + [30] + col_widths[4:];
            # For output to the file, we write the whole line each time
            # Add the initial entry fields here
            # This will also be used for some status updates where the whole message needs to be printed
            # to the screen
            
            out_line = [ spacedOut(str(out_line[i]), term_col_widths[i]) for i in range(len(out_line)) ];
            file_line = [ spacedOut(str(file_line[i]), col_widths[i]) for i in range(len(file_line)) ];
            # Compile both the truncated and the full status update

            if full_update:
                sys.stdout.write("".join(file_line) + "\n");
                sys.stdout.flush();
            else:         
                sys.stdout.write("\b" * 40);
                sys.stdout.write("".join(out_line) + "\n");
                sys.stdout.flush();
            # For full updates, print the full line to the screen
            # For others, delete the "In progress..." column and update the same status line
            
            printWrite(logfilename, 3, "".join(file_line));
            # Write the full line to the file.

        # The final status entry

        #####

    return cur_time;

#############################################################################
