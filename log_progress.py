import logging
import datetime

def setup_logger(log_file):
    # Initialize logging
    log_file = f"{log_file}.{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log.txt"
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s   %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

def log(message):
    logging.info(message)
    print(message)

if __name__ == "__main__":
    # Example of using the logger
    log_file = f"log_{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.txt"
    setup_logger(log_file)
    log("This is an example log message.")
