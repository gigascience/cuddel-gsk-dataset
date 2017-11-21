import os
from isatools import isatab

my_json_report = isatab.validate(open(os.path.join('./metadata-isa/', 'i_gsk_longitudinal.txt')))
