import numpy as np, random, operator, pandas as pd, matplotlib.pyplot as plt
from tools.geofunc import GeoFunc
from tools.show import PltFunc
from tools.nfp import NFP
from tools.data import getData
from tools.packing import PackingUtil,NFPAssistant,PolyListProcessor,Poly
from tools.bottom_left_fill import BottomLeftFill
import json
from shapely.geometry import Polygon,mapping
from shapely import affinity
import csv
import time
import multiprocessing
import datetime
import random
import copy

# fitness
def packingLength(poly_list,history_index_list,history_length_list,width,**kw):
    polys=PolyListProcessor.getPolysVertices(poly_list)
    index_list=PolyListProcessor.getPolyListIndex(poly_list)
    length=0
    check_index=PolyListProcessor.getIndex(index_list,history_index_list)
    first_call = False
    if 'is_first_call' in kw:
        first_call = kw['is_first_call']
    # if check_index>=0:
    #     length=history_length_list[check_index]
    # else:
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

class DE(object):

    def __init__(self, PopSize,Iteration,width,polys,all_rotation, **kw):
        self.width = width
        self.polys = polys
        self.nfp_assistant = None
        if 'nfp_asst' in kw:
            self.nfp_assistant = kw['nfp_asst']
        self.all_rotation=all_rotation  # 允许旋转的角度
        self.Iteration=Iteration  #  迭代次数
        self.PopSize=PopSize # 种群大小
        # self.Iteration=Iteration*PopSize # 最大函数评估次数
        self.Nmin = 4  # 最小种群大小
        self.H = PopSize  # 历史存档大小
        self.MCR = 0.5 * np.ones(self.H)  # 初始化交叉概率历史存档
        self.MF = 0.5 * np.ones(self.H)  # 初始化变异因子历史存档
        self.Fsigma = 0.1  # 变异因子的标准差
        N = PopSize  # 当前种群大小
        self.FES = N  # 当前函数评估次数
        self.area=0 # 计算所有多边形的面积
        for i in range(len(polys)):
            poly_shape = Polygon(polys[i])
            self.area += poly_shape.area # 面积
        # 初始化，要把每个cur_poly_list中的矩形的下标乱序一下
        self.cur_poly_list = [] 
        for i in range(self.PopSize): 
            # 旋转所有多边形
            for i in range(len(polys)):
                poly_shape = Polygon(polys[i])
                # 随机选择一个旋转角度
                rotation_angle = random.choice(all_rotation)
                rotated_poly = affinity.rotate(poly_shape, rotation_angle)
                coords = list(rotated_poly.exterior.coords)[:-1]
                polys[i] = [[x, y] for x, y in coords]
                random.shuffle(polys) # 乱序
            poly_list = PolyListProcessor.getPolyObjectList(polys, all_rotation)
            
            self.cur_poly_list.append(poly_list) # 当前的序列
        self.new_poly_list=self.cur_poly_list # 生成新的序列
        
        self.history_index_list=[] # 运行过的index序列
        self.history_length_list=[] # 运行结果
        
        self.NFPAssistant=NFPAssistant(PolyListProcessor.getPolysVertices(poly_list),get_all_nfp=True)

        self.run()
  
    def run(self):
        initial_length = [0] * self.PopSize
        # for length都计算出来 , self.cur_poly_list[i], 要用for i PopSize 取 i  
        for i in range(self.PopSize): 
            initial_length[i]=packingLength(self.cur_poly_list[i],self.history_index_list,self.history_length_list,self.width, is_first_call=True)

        global_lowest_length_list = [] # 记录每代最最小长度，理论上会下降
        temp_lowest_length_list= [] # 每代平衡长度

        global_lowest_length=min(initial_length) # 全局最小长度
        print(global_lowest_length)
        min_index = initial_length.index(global_lowest_length)
        global_best_list = copy.deepcopy(self.cur_poly_list[min_index]) # 用于记录历史上最好
        
        temp_best_list=copy.deepcopy(self.cur_poly_list[min_index]) # 局部最低排样
        temp_lowest_length=min(initial_length) # 局部最低长度

        q = 1  # 历史存档更新索引
        it=1
        # 开始循环寻找
        while self.Iteration>it:
            print("当前代数:",it)
            SCR = []  # 成功交叉概率
            SF = []  # 成功变异因子
            pmin = 2 / PopSize  # 最小选择概率
            cur_length_list=initial_length = [0] * self.PopSize
            old_lowest_length=global_lowest_length # 统计未更改次数
            for i in range(self.PopSize): 
                cur_length_list[i]=packingLength(self.cur_poly_list[i],self.history_index_list,self.history_length_list,self.width,NFPAssistant=self.NFPAssistant)
            cur_minlength=min(cur_length_list)
            # 种群循环
            for i in range(self.PopSize): 
                
                r = np.random.randint(self.H)  # 随机选择一个历史存档索引
                CR = self.MCR[r] + 0.1 * np.random.randn()  # 生成交叉概率
                CR = max(0, min(1, CR))  # 将交叉概率限制在 [0, 1] 范围内
                F = self.MF[r] + self.Fsigma * np.tan(np.pi * (np.random.rand() - 0.5))  
                F = min(1, F)  # 限制在 [0, 1] 范围内
                while F <= 0:  # 如果小于等于 0，重新生成
                    F = self.MF[r] + self.Fsigma * np.tan(np.pi * (np.random.rand() - 0.5))
                    F = min(1, F)
                
                #变异
                rs = list(range(len(self.cur_poly_list)))
                # 打乱列表顺序
                random.shuffle(rs)
                self.new_poly_list[i]=PolyListProcessor.mutate(temp_best_list,self.cur_poly_list[rs[0]],
                        self.cur_poly_list[rs[1]],self.cur_poly_list[rs[2]],self.cur_poly_list[rs[3]],F)
                
                # 交叉
                for j in range(int(len(global_best_list))):
                    if random.random()< CR:
                        self.new_poly_list[i]=PolyListProcessor.cross(self.new_poly_list[i],self.cur_poly_list[i],j,self.all_rotation)

                new_length=packingLength(self.new_poly_list[i],self.history_index_list,self.history_length_list,self.width,NFPAssistant=self.NFPAssistant)
                delta_length = new_length-cur_minlength

                if delta_length < 0: # 当代长度更小则接受
                    temp_best_list = self.cur_poly_list[i] = copy.deepcopy(self.new_poly_list[i]) 
                    temp_lowest_length=new_length # 修改为新的长度
                    cur_minlength=new_length

                    if new_length<global_lowest_length: # 如果新的长度小于最小的长度则修改最低长度
                        global_lowest_length=new_length
                        global_best_list=copy.deepcopy(self.new_poly_list[i])
                        cur_minlength=new_length
                        SCR.append(CR)  # 记录成功交叉概率
                        SF.append(F)  # 记录成功变异因子

                    elif np.random.random() < np.exp(-delta_length / self.Iteration): # 按照一定概率修改，并作为下一次检索基础
                        self.poly_list=copy.deepcopy(self.new_poly_list[i])
                        cur_minlength=new_length
                    else:
                        pass # 否则不进行修改
            print(global_lowest_length)
            # print("当前温度最低长度:",temp_lowest_length)
            # print("最低长度:",global_lowest_length)

            # 更新 MCR 和 MF
            if len(SCR) > 0 and len(SF) > 0:  # 如果有成功的交叉和变异因子
                b = np.abs(global_lowest_length - cur_length_list[q])  # 计算适应度差异
                sumf = np.sum(b)  # 计算适应度差异的总和
                if sumf > 0:  # 检查 sumf 是否为零
                    w = b / sumf  # 计算权重
                    self.MCR[q] = np.sum(w * np.array(SCR))  # 更新交叉概率历史存档
                    self.MF[q] = np.sum(w * np.array(SF)**2) / np.sum(w * np.array(SF))  # 更新变异因子历史存档
                    q += 1  # 更新历史存档索引
                    if q >= self.H:  # 如果索引超出范围，重置为 0
                        q = 0
            # self.cur_poly_list=copy.deepcopy(temp_best_list) # 某代检索结束后取该代最优值
            it+=1
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
        
        # 长度变化
        # self.showBestResult(temp_lowest_length_list,global_lowest_length_list)
    
    def showBestResult(self,list1,list2):
        plt.figure(1)
        plt.subplot(311)
        plt.plot(list1)#每代平衡路径长度
        plt.subplot(312)
        plt.plot(list2)#每代最好路径长度
        plt.grid()
        plt.show() 

if __name__=='__main__':
    starttime = datetime.datetime.now()
    #4:dighe1:100,6：fu：38，10：mao：510
    width=38
    polys = getData(6)
    # all_rotation = [0,90,180,270] # 旋转
    all_rotation = [0,180] # 旋转
    Iteration=300
    PopSize=30
    nfp_assistant=NFPAssistant(polys, store_nfp=False, get_all_nfp=True, load_history=True)

    DE(PopSize,Iteration,width,polys,all_rotation,nfp_asst=nfp_assistant)



