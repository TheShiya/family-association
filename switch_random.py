f = open("id_row_dup.txt", "r")
n = open("new_id_row.txt", "w")

for line in f:
	#print(type(line))
	if "0" in line:
		#print("no")
		n.write("1")
		n.write("\n")
	if "1" in line:
		n.write("0")
		n.write("\n")

n.close()
