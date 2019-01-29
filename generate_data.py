import numpy as np
nPart = 100

mass = (np.random.randn(nPart)+3)/3.0
x = np.random.randn(nPart)*10 +100
y = np.random.randn(nPart)*10 +100
z = np.random.randn(nPart)*10 +100
vx = np.random.randn(nPart)+10
vy = np.random.randn(nPart)+10
vz = np.random.randn(nPart)+10
softening = np.full(nPart,np.mean(x))
potential = np.zeros(nPart)

myFile = open("test_data.txt",'w+')
myFile.write("%d 0 %d\n" % (nPart,nPart))
myFile.write("".join(str(elem)+"\n" for elem in mass))
myFile.write("".join(str(elem)+"\n" for elem in x))
myFile.write("".join(str(elem)+"\n" for elem in y))
myFile.write("".join(str(elem)+"\n" for elem in z))
myFile.write("".join(str(elem)+"\n" for elem in vx))
myFile.write("".join(str(elem)+"\n" for elem in vy))
myFile.write("".join(str(elem)+"\n" for elem in vz))
myFile.write("".join(str(elem)+"\n" for elem in softening))
myFile.write("".join(str(elem)+"\n" for elem in potential))
myFile.close()
