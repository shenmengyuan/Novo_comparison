with open("groups.txt", 'r') as lines:
	for line in lines:
		if (("Aro" in line) & ("Len" in line) & ("Mbe" in line) & ("P6w" in line) & ("Pp1" in line) & \
			("Bar" in line) & ("Fuk" in line) & ("Lin" in line) & ("Pan" in line) & \
			("Pen" in line) & ("Res" in line) & ("Spb" in line) & ("thn" in line) & ("SCN66" in line) & \
	          	("SCN63" in line) & ("sp637" in line) & ("ST904" in line) & ("KN65" in line) & \
			("PC22D" in line) & ("Chol11" in line) & ("subter" in line) & ("napht" in line)):
			print line,
