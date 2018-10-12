## Validation of ISA-Tab metadata produced manually and available in folder 'metadata-isa'

import os
import logging
from isatools import isatab

my_json_report = isatab.validate(open(os.path.join('./metadata-isa/', 'i_gsk_longitudinal.txt')), log_level=logging.WARN)

