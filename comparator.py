import sys

def determine(lstELM) :
    if(lstELM == 'Remove' or  lstELM == 'Follow' or lstELM == 'Hashtags' or lstELM=='People' or lstELM=="·" or lstELM=='Following') :
        return True
    elif(lstELM.find(" ") != -1) :
        return True
    return False

names = ["boneS10", "Emilia_923", "ldoor", "af_5_k101", "Serena", "audikw_1"]
z=int(sys.argv[1])
file1 = open('/home/selin/3way-Par-Results/' + names[z] + '/result.txt')
data = file1.read().split('\t')
file2 = open('/home/selin/3way-Seq-Results/' + names[z] + '/result.txt')
data2 = file2.read().split('\t')
file1.close()
file2.close()
if(len(data) != len(data2))  :
    print("Size error.")

for i in range(len(data)):
    # matching line1 from both files
    try:
        if abs(float(data[i]) - float(data2[i])) > 0.01:
            print("\tLine-", i, end='')
            print("\tPar:", float(data[i]), end='')
            print("\tSeq:", float(data2[i]), end='')
            print("\n")

    except ValueError:
        print("\tLine-", i, end='')
        print(data[i])
        print(data2[i])
