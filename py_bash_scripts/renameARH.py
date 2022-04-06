import os
import sys

"""takes all t_TIME_pRANK.arh files and saves them in a folder 
corresponding to different times while renaming them to pRANK.arh"""

def main():
    filepath = "/scratch/ws/1/haja565a-workspace1/collision/out151/"
    filepath = sys.argv[1]
    filepath += "/"

    prefix = "t_"

    for filename in os.listdir(filepath):
        if(prefix == filename[:2]):
            #print(filename)
            src =filepath+filename
            new_name = filename[9:];
            #print(new_name)
            newpath = filepath + "time" + filename[1:8]+ "/"
            if not os.path.exists(newpath):
                os.mkdir(newpath)
            dst = newpath + new_name
            print("***")
            print("renaming:" + src)
            print("to: " + dst)
            os.rename(src, dst)
    

if __name__ == '__main__':
    main()