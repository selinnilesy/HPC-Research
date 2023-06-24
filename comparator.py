import sys

names = ["boneS10", "Emilia_923", "ldoor", "af_5_k101", "Serena", "audikw_1"]
z=int(sys.argv[1])
file1 = open('/home/selin/3way-Par-Results/' + names[z] + '/result.txt')
data = file1.read().split('\t')
file2 = open('/home/selin/3way-Seq-Results/' + names[z] + '/result.txt')
#file2 = open('/home/selin/Seq-Results/' + names[z] + '/banded/result.txt')
data2 = file2.read().split('\t')
file1.close()
file2.close()
if(len(data) != len(data2))  :
    print("Size error.")
else :
    print("Output elements size:",len(data))

for i in range(len(data)):
    # matching lines from both files
    try:
        if abs(float(data[i]) - float(data2[i])) > 0.000001:
            print("\tLine-", i)
            print("\tPar:", float(data[i]))
            print("\tSeq:", float(data2[i]))
            print("\n")
            break

    except ValueError:
        print("\tLine-", i)
        print(data[i])
        print(data2[i])
