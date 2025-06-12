from tools.geofunc import GeoFunc
import pandas as pd
import json

def getData(index): #0-12,13-19
    name=["ga","albano","blaz1","blaz2","dighe1","dighe2","fu","han","jakobs1","jakobs2","mao","marques","shapes",
          "shirts","swim","trousers","1-cavity","1-degrade","2-degrade","1-shirt","shirts_11","5-shirt","trousers20"]
    print("开始处理",name[index],"数据集")
    scale=[1,1,1,1,1,1,1,1,1,1,0.2,1,1,1,1,1,1,1,1,1,1,1,1]
    # print("缩放",scale[index],"倍")
    df = pd.read_csv("data/"+name[index]+".csv")
    polygons=[]
    for i in range(0,df.shape[0]):
        for j in range(0,df['num'][i]):
            poly=json.loads(df['polygon'][i])
            GeoFunc.normData(poly,scale[index])
            polygons.append(poly)
    return polygons
