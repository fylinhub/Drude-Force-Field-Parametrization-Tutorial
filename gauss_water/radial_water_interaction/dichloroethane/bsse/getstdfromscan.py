#!/usr/bin/python

from sys import argv

entrada = open(argv[1],'r')

EOF = 0
cadena = ''
n = 0
while (cadena.find('Normal termination') == -1): 
	while (cadena.find('Standard orientation') == -1):
		cadena = entrada.readline()
	for i in range(5):
		cadena = entrada.readline()
	lista = []
	nombre = 'geom-' + str(n) + '.xyz'
	while (cadena.find('---') == -1):
#		print cadena.split()[1],cadena.split()[3],cadena.split()[4],cadena.split()[5]
		lista.append([cadena.split()[1],cadena.split()[3],cadena.split()[4],cadena.split()[5]])
		cadena = entrada.readline()
	salida = open(nombre,'w')
	salida.write(str(len(lista)))
	salida.write('\n')
	salida.write('\n')
	for i in range(len(lista)):
		linea = '%s %s %s %s\n' % (lista[i][0],lista[i][1],lista[i][2],lista[i][3])
		#print linea
		salida.write(linea)
	salida.close()
	n = n + 1
entrada.close()
