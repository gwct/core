#############################################################################
# Python logging configuration
#############################################################################

import sys, os
import logging
import __main__

#############################################################################

# Example of your ColoredFormatter
class ColoredFormatter(logging.Formatter):
    RESET = "\033[0m"
    COLOR_MAP = {
        logging.DEBUG: "\033[35m",      # Purple (Magenta) for DEBUG
        logging.INFO: "\033[36m",       # Cyan for INFO
        logging.WARNING: "\033[33m",    # Yellow for WARNING
        logging.ERROR: "\033[31m",      # Red for ERROR
        logging.CRITICAL: "\033[1;31m"    # Bold red for CRITICAL
    }

    def format(self, record):
        color = self.COLOR_MAP.get(record.levelno, self.RESET)
        message = super().format(record)
        return f"{color}{message}{self.RESET}"

# Define the filter that will block records flagged as file_only.
def no_file_only(record):
    return not getattr(record, 'file_only', False)

def configureLogging(log_level: str = "INFO", log_verbosity: str = "BOTH", log_filename: str = "", logger_name: str = "core_logger", 
                        show_name: bool = False, overwrite_log_file: bool = False, usage: bool = False) -> None:
    """
    configureLogging(log_level='INFO', log_verbosity='BOTH', log_filename='', logger_name='core_logger', show_name=False, usage=False)

    Setup a logger with customizable output options including file, terminal (stdout/stderr), and log formatting.
    Logging colors determined by the ColoredFormatter class.

    Parameters:
    ----------
    log_level : str, optional
        The logging level. One of: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'. Default: 'INFO'.
    
    log_verbosity : str, optional
        Determines where logs are output. One of: 'BOTH', 'FILE', 'SCREEN'. Default: 'BOTH'.
    
    log_filename : str, optional
        The file to write logs to. If empty and verbosity includes 'FILE', a default filename is used (based on script name). Default: ''.
    
    logger_name : str, optional
        Name of the logger instance to configure. Default: 'core_logger'.
    
    show_name : bool, optional
        Whether to include logger name in the log output format. Default: False.

    overwrite_log_file : bool, optional
        If True, overwrites the log file if it already exists, otherwise appends to it. Default: False.
    
    usage : bool, optional
        If True, prints this usage message and exits without configuring logging.

    Returns:
    -------
    str
        The name of the logger configured (useful for referencing later).

    Example Usage:
    --------------
    import logging
    logger_name = configureLogging(log_level='INFO', log_verbosity='BOTH')
    my_logger = logging.getLogger(logger_name)
    my_logger.error("An exception occurred", exc_info=True)
    my_logger.info("This only goes to the file", extra={"file_only": True})

    """

    if usage:
        print(configureLogging.__doc__);
        sys.exit("Exiting.");

    log_verbosity = log_verbosity.upper();
    if log_verbosity not in ["BOTH", "FILE", "SCREEN"]:
        raise ValueError(f"Error in loginit: Invalid log verbosity: {log_verbosity}. Choose from BOTH, FILE, SCREEN.")
    # Some error checking

    core_logger = logging.getLogger(logger_name);
    # Set the logger name to the specified name

    if log_filename == "" and log_verbosity in ["FILE", "BOTH"]:
        log_filename = os.path.splitext(os.path.basename(__main__.__file__))[0] + ".log"
    # Default log filename if not provided and verbosity is FILE or BOTH

    if show_name:
        log_format = '[ %(asctime)s  - %(name)s - %(levelname)s ] %(message)s';
    else:
        log_format = '[ %(asctime)s - %(levelname)s ] %(message)s';
    date_format = '%Y-%m-%d %H:%M:%S'
    # Format for the log messages
    
    # Check if logger has handlers already to avoid duplicate handlers
    if not core_logger.hasHandlers():
        if log_verbosity in ["BOTH", "FILE"]:
            file_mode = 'w' if overwrite_log_file else 'a';
            handler_file = logging.FileHandler(log_filename, mode=file_mode);
            handler_file.setFormatter(logging.Formatter(fmt=log_format, datefmt=date_format))
            core_logger.addHandler(handler_file)
        # Add file handler if specified

        if log_verbosity in ["BOTH", "SCREEN"]:
            # Create a handler for DEBUG and INFO messages that prints to stdout
            handler_stdout = logging.StreamHandler(sys.stdout)
            handler_stdout.setFormatter(ColoredFormatter(fmt=log_format, datefmt=date_format))
            # Only allow messages with level less than WARNING (DEBUG, INFO)
            handler_stdout.addFilter(lambda record: record.levelno < logging.WARNING)
            # And add our extra filter to skip file-only records
            handler_stdout.addFilter(no_file_only)
            core_logger.addHandler(handler_stdout)

            # Create a separate handler for WARNING and above that prints to stderr
            handler_stderr = logging.StreamHandler(sys.stderr)
            handler_stderr.setFormatter(ColoredFormatter(fmt=log_format, datefmt=date_format))
            # Only allow messages with level WARNING or higher
            handler_stderr.addFilter(lambda record: record.levelno >= logging.WARNING)
            # Also add our filter here to skip file-only records
            handler_stderr.addFilter(no_file_only)
            core_logger.addHandler(handler_stderr)
    
    # Set the logging level
    log_level = log_level.upper()
    level_map = {
        'CRITICAL': logging.CRITICAL,
        'ERROR': logging.ERROR,
        'WARNING': logging.WARNING,
        'INFO': logging.INFO,
        'DEBUG': logging.DEBUG,
        'NOTSET': logging.NOTSET
    }
    try:
        core_logger.setLevel(level_map[log_level])
        core_logger.debug(f"Logging level set to {log_level}")
    except KeyError:
        raise ValueError(f"Error in loginit: Invalid logging level: {log_level}. Choose from {list(level_map.keys())}")

    return logger_name;
    # Return the logger name if the default is used

#############################################################################

if __name__ == "__main__":
    # Configure the logger
    logger_name = configureLogging("debug", "both", usage=True)
    my_logger = logging.getLogger(logger_name)

    # Test logging at different levels
    my_logger.debug("This is a debug message.")
    my_logger.info("This is an info message.")
    my_logger.warning("This is a warning message.")
    my_logger.error("This is an error message.")
    my_logger.critical("This is a critical message.")

    # Test logging with variables
    test_variable = "test_value"
    my_logger.info(f"Logging a variable: {test_variable}")

    # Test logging exceptions
    try:
        1 / 0
    except ZeroDivisionError as e:
        my_logger.error("An exception occurred", exc_info=True)
