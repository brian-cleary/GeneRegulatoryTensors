import sys,getopt,os


if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'p:n:',["outfile="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-p'):
			prefix = arg
		elif opt in ('-n'):
			n = int(arg)
	x = 0
	f = open(prefix + '.matches.txt')
	g = open('%s.%d-%d.matches.txt' % (prefix,x,x+n-1),'w')
	for line in f:
		ls = line.strip().split()
		if int(ls[2]) >= x+n:
			g.close()
			os.system('gzip %s.%d-%d.matches.txt' % (prefix,x,x+n-1))
			x += n
			g = open('%s.%d-%d.matches.txt' % (prefix,x,x+n-1),'w')
		g.write(line)
	f.close()
	g.close()
	os.system('gzip %s.%d-%d.matches.txt' % (prefix,x,x+n-1))
	os.system('rm %s.matches.txt' % prefix)