#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK

"""
PDF grepping tool
Search for string pattern in given PDF if contains text.
Use pypdfocr to convert image pdf to text-based without altering layout
Example usage:
    pdfgrep.py -irl Bank *.pdf
"""


from __future__ import print_function

__author__ = "Amro Diab"
__author_email__ = "adiab@linuxmail.org"

import argparse
import os
import re
import sys

from glob import glob
import argcomplete
import magic
import pdftotext

APP_NAME = sys.argv[0].split('/')[-1]
def pdfgrep(grep, filenames, ignoreCase=False, recursive=False, listFiles=False, color=False, num=False):
	"""
	Main function, fetch command line arguments and call do_grep with given
	attibutes
	"""

	if not filenames:
		print("No filename provided, pdfgrep terminating.")
		exit(1)
		#filenames = sorted(glob('*'))
		
	for filename in filenames:
		for output in \
			do_grep(filename, grep, ignoreCase=ignoreCase,
					num=num, listFiles=listFiles,
					recursive=recursive, color=color):
			print(output)

def do_grep(filename, grep, **kwargs):
	""" Perform the pattern matching in given files """

	for arg in ['ignoreCase', 'listFiles', 'num', 'recursive', 'color']:
		kwargs.setdefault(arg, False)

	# set ignore case flag
	re.IGNORECASE = 0 if not kwargs['ignoreCase'] else 2

	# recurse through directory structure if filename is a dir
	if os.path.isdir(filename) and kwargs['recursive'] is True:
		for r_file in glob(filename + '/*'):
			for output in \
				do_grep(r_file, grep, ignoreCase=kwargs['ignoreCase'],
						listFiles=kwargs['listFiles'],
						num=kwargs['num'], recursive=kwargs['recursive'],
						color=kwargs['color']):

				print(output)
		return

	elif os.path.isdir(filename) and kwargs['recursive'] is False:
		# ignore directory and return if recursive flag is not set
		return
	elif not os.path.isfile(filename):
		sys.stderr.write("{0}: {1}: No such file or directory\n"
						 .format(APP_NAME, filename))
		return

	elif magic.from_file(filename, mime=True) != 'application/pdf':
		sys.stderr.write("{0}: Not a pdf file: {1}\n "
						 .format(APP_NAME, filename))
		return

	try:
		pdf_file = open(filename, 'rb')
	except IOError:
		sys.stderr.write("{0}: {1}: No such file or directory\n"
						 .format(APP_NAME, filename))
		return

	try:
		read_pdf = pdftotext.PDF(pdf_file)

		# This will happen if file is malformed, or not a PDF
	except (pdftotext.Error, IOError):
		sys.stderr.write("{0}: Unable to read file: {1}\n "
						 .format(APP_NAME, filename))
		return
	except KeyboardInterrupt:
		return

	for page_num in range(0, len(read_pdf)):
		# attempt to read pages and split lines approprietly
		page = read_pdf[page_num]
		page_content = page.split('\n')

		# iterate through pages
		for line_num, line in enumerate(page_content):
			line = str(line.encode('utf-8'))
			if re.search(grep, line, re.IGNORECASE):
				if kwargs['listFiles']:
					# only print file names - return after first match found
					yield filename
					return
				# add page and line numbers if num flag set
				beg = ("page:{0}, line:{1}"
					   .format(page_num + 1, line_num + 1) if kwargs['num'] else "")

				if kwargs['color']:
					# highlight matching text in red
					red = '\033[31m'
					end = '\033[0m'
					text = re.compile(re.escape(grep), re.IGNORECASE)
					line = text.sub(red + grep + end, line)

				yield "{0}: {1} {2}".format(filename, beg, line.strip())
