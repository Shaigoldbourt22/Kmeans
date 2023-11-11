import sys
import numpy as np
import pandas as pd
import spkmeansmodule


def choose_centroids(vec_points, k):
    np.random.seed(0)
    index_arr = [0]*k
    n=len(vec_points)

    vec_cent=[]
    vec_d=[0]*n
    vec_p=[1/n]*len(vec_points)
    rand_num=np.random.choice(np.arange(0, n))
    index_arr[0]=rand_num
    cent=vec_points[rand_num]
    vec_cent.append(cent)
    for t in range(1,k):
        for i in range(n):
            vec_diff = [oc_norm(vec_points[i], vec_cent[j]) for j in range(t)]
            vec_d[i]=min(vec_diff)
        for z in range(n):
            vec_p[z]=vec_d[z]/(sum(vec_d))

        rand_num = np.random.choice(n , p=np.array(vec_p))
        cent=vec_points[rand_num]
        index_arr[t]=rand_num
        vec_cent.append(cent)
    str = ""
    for i in range(len(index_arr)):
        str += f"{index_arr[i]}"
        if i!= len(index_arr)-1:
            str+= ","
    print(str)
    return np.array(vec_cent)





def create_txt(txt_1):
    text1 = np.genfromtxt(txt_1, delimiter=",")
    text1 = pd.DataFrame(data=text1[0:,1:],index=text1[0:,0])
    return text_1



def oc_norm(point1, point2):
    sum_coor = 0
    for t in range(0, len(point1)):
        sum_coor = sum_coor+(point1[t]-point2[t])**2

    return sum_coor

def is_not_num(x):
    for i in range(len(x)):
        if not x[i].isdigit():
            return 1


def aplly_algo(points, k, goal):
    goal_dict = {"spk": 1, "wam" : 2 , "ddg" : 3, "lnorm" : 4, "jacobi" : 5}
    n = points.shape[0]

    if goal_dict[goal]==1 and k>0:
        centroids = choose_centroids(points,k)
    else:
        centroids = np.zeros((n,n))
    k_is_zero = (k==0)
    t_matrix = spkmeansmodule.fit(points.tolist(), centroids.tolist(), goal_dict[goal], k_is_zero)

    matrix = np.array(t_matrix)


    if goal_dict[goal]==1 and k == 0:
        new_k = matrix.shape[1]
        if new_k != 0 and k == 0:
            aplly_algo(matrix, new_k, goal)
        if k==0 and new_k==0:
            assert (1)

# main
array = sys.argv

if len(array) != 4:
    print("Invalid Input!")
    exit(1)


if is_not_num(array[1])==1:
    print("Invalid Input!")
    exit(1)



if int(array[1])<0:
    print("Invalid Input!")
    exit(1)

k=int(array[1])
goal=array[2]
file_name=array[3]
points = np.genfromtxt(file_name, delimiter=',')
if(k > points.shape[0] or k<0 or k == 1):
    print("Invalid Input!")
    exit(1)



aplly_algo(points,k, goal)


