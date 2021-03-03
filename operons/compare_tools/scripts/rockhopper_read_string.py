import sys

files = sorted(sys.argv[1:])

print( ",".join([files[ix*2]+ "%" + files[ix*2 + 1] for ix in range(len(files)//2)]) )