import csv

#Path to the input and ouput files
filename = "/home/margaret/mastersProject/vcf_files/new_file.csv"
outfile = "/home/margaret/mastersProject/vcf_files/fasta_file"

#Initializing empty lists
rows = []
fasta = []

#Opening the files for reading
with open(filename, 'r') as csvfile: 
    csvreader = csv.reader(csvfile)
    with open(outfile, 'a') as out:
    
        for row in csvreader:         #reads the csv file row after row and appends them to the initially empty list called rows
            rows.append(row)
            
        for row in rows[4:]:          #appends the 4th to last row to the empty list called fasta
            fasta.append(row)
            
        for list in fasta:            #opens the list containing the rows representing the isolates   
            for word in list:       #accesses each of the nested list(each isolate)
                if len(word) > 1:
                    header = '\n' + ">" + word + '\n'
                    out.write(word)              #prints the headers into the output file

                else:
                    if len(word) == 1:        #prints the sequences to the output file
                        out.write(word)
                