
import sys
import logging
logging.basicConfig(format='%(levelname)s:%(message)s')#, level=logging.DEBUG)


delay = float(sys.argv[1]) # seconds
type = sys.argv[2] #'event' or 'trajectory'
if type == 'event': spiking = True
elif type == 'trajectory': spiking = False
else: raise ValueError(f'option "{type}" is not valid; only "event" or "trajectory" are')
readname = sys.argv[3] #'fifa'
writname = sys.argv[4] #'fifb'
# TODO the roundoff can add a jitter of (sim dt) as it stands, make it work assuming dt is an integer timestep of usec...
delay = 0.010 

# s = 'StartA\n'
logging.debug(f'ReadyA {type} {readname} {writname}')
# buffering is 0 to overcome 'not seekable' and binary is set bc can't have unbuffered *text* IO
# and write is not set because then the pipe would remain open (for this process that asked) forever
# Open the pipe for reading BEFORE the one for writing, because the latter can't be opened until the file has been opened for reading on the other side, likewise for the other process...
# or use async IO TODO
with open(readname, 'rb',buffering=0) as fi:
	logging.debug('ReadyA')
	with open(writname, 'wb',buffering=0) as fo:
		logging.debug('ReadyA')
		i = 0
		if spiking:
			# logging.debug('waaaa!!')
			fo.write(f'{delay}\n'.encode("utf-8"))
		while fi:
			s = fi.readline()
			if not s:
				break
			# print(s)
			if not spiking and i == 0: 
				fo.write(s) # we still need that first line for eden to begin, parrot the first one that eden does emit. TODO a more elegant solution :/
			s = s.decode("utf-8")
			timestamp, sep, remainder = s.strip().partition(' ')
			# when there is only the timestamp+'\n' it is not split... strip it and add \n
			# or lstrip() + re.findall(r'^\S', st)[0] to preserve trailnig whitespace
			upds = str(float(timestamp)+delay)+sep+remainder+'\n'#+f'A {i}! \n'
			# print('upd', upds.encode("utf-8"))
			# logging.debug(f'A {i}')
			fo.write(upds.encode("utf-8"))
			i += 1
