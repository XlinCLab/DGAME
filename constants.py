import logging
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(levelname)s: %(message)s')

AUDIO_FILE_PATTERN = re.compile(r"(\d+)_words_(\d+)\.csv")
