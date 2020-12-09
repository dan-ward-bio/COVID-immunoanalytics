#! /usr/bin/env python3

# First argument is window size
# Second arugment is string

import sys
windowsize = sys.argv[1] 
intwindowsize = int(windowsize)
String = str(sys.argv[2]) 
FragmentsOfSize = [String[i:i+intwindowsize] for i in range(len(String))]
for fragment in FragmentsOfSize:
	print(fragment)