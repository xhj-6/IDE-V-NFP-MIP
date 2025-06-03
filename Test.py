import numpy as np, random, operator, pandas as pd, matplotlib.pyplot as plt
from tools.geofunc import GeoFunc
from tools.show import PltFunc
from tools.nfp import NFP
from tools.data import getData
from tools.packing import PackingUtil,NFPAssistant,PolyListProcessor,Poly
from heuristic import TOPOS,BottomLeftFill
import json
from shapely.geometry import Polygon,mapping
from shapely import affinity
import csv
import time
import multiprocessing
import datetime
import random
import copy

def packingLength(poly_list,history_index_list,history_length_list,width,**kw):
    polys=PolyListProcessor.getPolysVertices(poly_list)
    index_list=PolyListProcessor.getPolyListIndex(poly_list)
    length=0
    check_index=PolyListProcessor.getIndex(index_list,history_index_list)
    first_call = False
    if 'is_first_call' in kw:
        first_call = kw['is_first_call']
    if check_index>=0:
        length=history_length_list[check_index]
    else:
        try:
            # if 'NFPAssistant' in kw:
            #     blf=BottomLeftFill(width,polys,NFPAssistant=kw['NFPAssistant'])
            #     # blf.showAll()
            #     length=blf.contain_length
            # else:
            length=BottomLeftFill(width,polys, is_first_call=first_call).contain_length
        except:
            # print('出现Self-intersection')
            length=99999
        
        history_index_list.append(index_list)
        history_length_list.append(length)
    return length

class SA(object):
    '''
    Simulating Annealing + Bottom Left Fill
    Reference:....
    '''
    def __init__(self, Iteration, width, polys,all_rotation, **kw):
        self.Iteration=Iteration  #  迭代次数
        self.width = width
        self.polys = polys
        self.nfp_assistant = None
        if 'nfp_asst' in kw:
            self.nfp_assistant = kw['nfp_asst']
        self.all_rotation=all_rotation  # 允许旋转的角度
        self.loop_times=1 # 内循环次数

        self.area=0
        # 旋转所有多边形
        for i in range(len(polys)):
            poly_shape = Polygon(polys[i])
            self.area += poly_shape.area # 面积
            # 随机选择一个旋转角度
            rotation_angle = random.choice(self.all_rotation)
            rotated_poly = affinity.rotate(poly_shape, rotation_angle)
            coords = list(rotated_poly.exterior.coords)[:-1]
            polys[i] = [[x, y] for x, y in coords]

        # 打乱polys的顺序
        random.shuffle(polys)

        poly_list = PolyListProcessor.getPolyObjectList(polys, self.all_rotation)
        
        self.cur_poly_list=poly_list # 当前的序列
        self.new_poly_list=poly_list # 生成新的序列

        self.history_index_list=[] # 运行过的index序列
        self.history_length_list=[] # 运行结果
        
        self.NFPAssistant=NFPAssistant(PolyListProcessor.getPolysVertices(poly_list),get_all_nfp=True)

        self.run()
    
    def newPolyList(self):
        choose_id = int(random.random() * len(self.new_poly_list))
        '''进行交换和旋转的操作.'''
        if random.random()<=0.6:
            self.new_poly_list=PolyListProcessor.randomSwap(self.cur_poly_list,choose_id)
        else:
            self.new_poly_list=PolyListProcessor.randomRotate(self.cur_poly_list,self.all_rotation,choose_id)
  
    def run(self):
        # initial_length=PolyListProcessor.getLength(self.width,self.area,self.cur_poly_list)
        initial_length=packingLength(self.cur_poly_list,self.history_index_list,self.history_length_list,self.width, is_first_call=True)
        print(initial_length)
        global_lowest_length_list = [] # 记录每代最最小长度，理论上会下降
        temp_lowest_length_list= [] # 每代平衡长度

        global_best_list = copy.deepcopy(self.cur_poly_list) # 用于记录历史上最好蓄力
        global_lowest_length=initial_length # 全局最低高度
        
        temp_best_list=copy.deepcopy(self.cur_poly_list) # 局部最低排样
        temp_lowest_length=initial_length # 局部最低长度

        it=1
        # 开始循环寻找
        while self.Iteration>it:
            print("当前代数：",it)
            old_lowest_length=global_lowest_length # 统计未更改次数
            # cur_length=PolyListProcessor.getLength(self.width,self.area,self.cur_poly_list)
            cur_length=packingLength(self.cur_poly_list,self.history_index_list,self.history_length_list,self.width,NFPAssistant=self.NFPAssistant)

            # 进行一定次数的寻找 
            for i in range(self.loop_times): 
                self.newPolyList()
                # new_length=PolyListProcessor.getLength(self.width,self.area,self.new_poly_list)
                new_length=packingLength(self.new_poly_list,self.history_index_list,self.history_length_list,self.width,NFPAssistant=self.NFPAssistant)
                delta_length = new_length-cur_length

                if delta_length < 0: # 当前如果高度更低则接受
                    temp_best_list = self.cur_poly_list = copy.deepcopy(self.new_poly_list) 
                    temp_lowest_length=new_length # 修改为新的高度
                    cur_length=new_length

                    if new_length<global_lowest_length: # 如果新的高度小于最低的高度则修改最低高度
                        global_lowest_length=new_length
                        global_best_list=copy.deepcopy(self.new_poly_list)
                        cur_length=new_length

                    elif np.random.random() < np.exp(-delta_length / self.temp_now): # 按照一定概率修改，并作为下一次检索基础
                        self.poly_list=copy.deepcopy(self.new_poly_list)
                        cur_length=new_length
                    else:
                        pass # 否则不进行修改
                it=it+1
            # 每代的结果
            print("当前最低长度:",temp_lowest_length)
            print("最低长度:",global_lowest_length)

            # 迭代5次后输出结果
            # polys=PolyListProcessor.getPolysVertices(global_best_list)
            # BottomLeftFill(width,polys).showAll()

            self.cur_poly_list=copy.deepcopy(temp_best_list) # 某代检索结束后取该代最优值
            
            global_lowest_length_list.append(global_lowest_length) # 全局的在每代的最低长度，理论上一直在降低
            temp_lowest_length_list.append(temp_lowest_length) # 每代的最低长度
            
            
        # print('局部最优的序列:',temp_best_list)
        print('局部最优长度:',temp_lowest_length)
        
        print('最好序列长度:',global_lowest_length)
        
        endtime = datetime.datetime.now()
        print (endtime - starttime)

        # print('best_list:')
        # best_list=[]
        # for p in global_best_list: 
        #     best_list.append(p.poly)
        # print(best_list)
        print(PolyListProcessor.getLength1(self.width,self.area,global_best_list))
        # PolyListProcessor.showPolyList(self.width,self.area,global_best_list)
        
        # 长度变化
        # self.showBestResult(temp_lowest_length_list,global_lowest_length_list)
    
    def showBestResult(self,list1,list2):
        plt.figure(1)
        # plt.subplot(311)
        # plt.plot(list1)#每代平衡路径长度
        # plt.subplot(312)
        plt.plot(list2)#每代最好路径长度
        # plt.grid()
        plt.show() 

if __name__=='__main__':
    starttime = datetime.datetime.now()
    #4:dighe1:100,6：fu：38，10：mao：510, 16:1-cavity:20, 17:1-degrade:20 18:2-degrade:15，19：1-shirt,20:shirts_11:80
    Iteration=2
    width=80
    polys = getData(20)
    # all_rotation = [0,90,180,270] # 旋转
    all_rotation = [0,180] # 旋转
    # all_rotation = [0] # 旋转
    # 旋转所有多边形
    for i in range(len(polys)):
        poly_shape = Polygon(polys[i])
        # 随机选择一个旋转角度
        rotation_angle = random.choice(all_rotation)
        rotated_poly = affinity.rotate(poly_shape, rotation_angle)
        coords = list(rotated_poly.exterior.coords)[:-1]
        polys[i] = [[x, y] for x, y in coords]

    # 打乱polys的顺序
    # random.shuffle(polys)

    poly_list = PolyListProcessor.getPolyObjectList(polys, all_rotation)
    
    nfp_assistant=NFPAssistant(polys, store_nfp=False, get_all_nfp=True, load_history=True)

    SA(Iteration,width,polys,all_rotation,nfp_asst=nfp_assistant)

    # endtime = datetime.datetime.now()
    # print (endtime - starttime)
