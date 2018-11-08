#!/usr/bin/env python


import sys	 
import traceback 


###================ Galaxyy mods below	v34a v0 ============
def main():


	# check the command line
	sys.stdout.write(' starting set invalid values to mean value. Arguments=')
	sys.stdout.write(str(sys.argv[1:])+'\n')

	try:
		inFile				= sys.argv[1]
		outFile				= sys.argv[2]
	
		colDelimiter  = '\t'
		
		Matrix = []
		fin = open( inFile, 'rU')
		indata = fin.read()
	        # split the file into lines
		#print 'fin', indata[1:]
		a = indata[:].split('\n')	
		cnt = -1
		for i in a:	 # for each row from matrix
			tmp = i.replace('\n','') ##v22 fix /n
			b = tmp.split(colDelimiter)
			tmp2 = []
			for j in b:
				tmp2.append(j)
			Matrix.append(tmp2)
	
	#find text in numeric cells and convert to 
		badvalue = 'false'
			#use header row for correct number of rows
		numrows = len(Matrix)
                print(numrows)
		numcols = len(Matrix[0])
                print(numcols)
		for rows in range(1, numrows ):
			for cols in range(1,numcols):
				if numcols <= len(Matrix[rows] ):
					try:
                                                Matrix[rows][cols] = float(Matrix[rows][cols]) # if non real number will cause error
					except:
						#sys.stdout.write('Illegal Value '+str(Matrix[i][cols])+ ' in row and column '+rows +cols +'\n' )
						print('Illegal value at row, column - ', rows+1,cols+1 )
						Matrix[rows][cols] = 99999   # temporarily set to number then later to the mean value
						#badvalue = 'true'
	 			
 				else:  
 					if rows < numrows -1: print('skip row', rows+1, 'not enough columns')
 					#if i == numrows-1: numrows = numrows -1

 	
 	
	#if non valid values above then change them to the COLUMN Mean value 
		if badvalue == 'false' and numrows >1 and numcols > 1: 
			for rows in range(1, numrows ):
# 				skip = 0
 				top  = 0
				if len(Matrix[rows]) > 1:
					for cols in range(1,numcols):
						try:
							top =top + float(Matrix[rows][cols])
						except:
							junk = 0
	 					
					meanN  = top / (numcols -1)
					
	 				for cols in range(1,numcols):
						try:
							Matrix[rows][cols] = str('%.3f' % (abs(float(Matrix[rows][cols]) -meanN)))
							#temp = Matrix[k][i]
						except:
							junk = 0
	 					
					#print('For row', rows+1,' the mean =', meanN )

		else:
			print( 'bad value or insufficient columns or rows ' , badvalue, numcols, numrows)
				
			# have mean so go back and set bad values to mean
		fout = open( outFile, 'w') 
                fout.write(str(str("Sample_Average"))+'\n')	
		for rows in range(0, numrows ):
			if len(Matrix[rows])> 2:
				try:
					for cols in range(numcols, numcols-1):
						fout.write( str(str(Matrix[rows][cols])) +'\t' )
					fout.write( str(str(Matrix[rows][cols]))+'\n' )
				except: junk = 0
			#if (rows < numrows -1 and len(Matrix[cols]) > 1): fout.write( '\n')
			#fout.write( '\n')
									
		
		fout.close()
		fin.close()
	except Exception as e:
		 print( 'Usage: python MeanCenter  failed', e )
		 print traceback.format_exc()
		 sys.exit(-1)



##		  
	print 'Success ', numrows, ' rows mean centered \n'
##
	return
##

if __name__ == "__main__": main()
