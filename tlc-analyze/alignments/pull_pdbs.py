import urllib2

durl = 'http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=\
pdb&compression=NO&structureId={0}'

pdbs = ['3RSW', '2HMB', '1HMR', '1HMS', '1HMT', '3NPO', '3NQ3', '3NQ9', '4DQ3',
        '3UEX', '3UEV', '3UEU', '3UEW', '2F73', '3STN', '3STK', '3STM', '3RZY',
        '3HK1', '3P6C', '3P6D', '3P6E', '3P6F', '3P6G', '3P6H', '2HZR', '2HZQ', 
        '3DTQ', '3DSZ', '2RCQ', '2RCT']

for name in pdbs:
    grab = durl.format(name)
    with open('{0}.pdb'.format(name), 'w') as pdb_file:
        pdb_file.write(urllib2.urlopen(grab).read())
