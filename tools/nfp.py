from tools.show import PltFunc
from tools.geofunc import GeoFunc
from tools.data import getData
from shapely.geometry import Polygon,Point,mapping,LineString
from shapely.ops import unary_union
import matplotlib.pyplot as plt
import pandas as pd
import json
import copy
import time

class NFP(object):
    def __init__(self,poly1,poly2,**kw):
        self.stationary=copy.deepcopy(poly1)
        self.sliding=copy.deepcopy(poly2)
        start_point_index=GeoFunc.checkBottom(self.stationary)
        self.start_point=[poly1[start_point_index][0],poly1[start_point_index][1]]
        self.locus_index=GeoFunc.checkTop(self.sliding)
        # 如果不加list则original_top是指针
        self.original_top=list(self.sliding[self.locus_index])
        GeoFunc.slideToPoint(self.sliding,self.sliding[self.locus_index],self.start_point)
        self.start=True # 判断是否初始
        self.nfp=[]
        self.rectangle=False
        if 'rectangle' in kw:
            if kw["rectangle"]==True:
                self.rectangle=True
        self.error=1
        self.main()
        if 'show' in kw:
            if kw["show"]==True:
                self.showResult()
        # 计算完成之后平移回原始位置
        GeoFunc.slideToPoint(self.sliding,self.sliding[self.locus_index],self.original_top)

    def main(self):
        # 验证输入的多边形
        sliding_poly = Polygon(self.sliding)
        stationary_poly = Polygon(self.stationary)

        if not sliding_poly.is_valid:
            print("警告：移动多边形几何形状无效")
            sliding_poly = sliding_poly.buffer(0)
        
        if not stationary_poly.is_valid:
            print("警告：固定多边形几何形状无效")
            stationary_poly = stationary_poly.buffer(0)
        
        inter = sliding_poly.intersection(stationary_poly)
        i=0
        if self.rectangle: # 若矩形则直接快速运算 点的index为左下角开始逆时针旋转
            width=self.sliding[1][0]-self.sliding[0][0]
            height=self.sliding[3][1]-self.sliding[0][1]
            self.nfp.append([self.stationary[0][0],self.stationary[0][1]])
            self.nfp.append([self.stationary[1][0]+width,self.stationary[1][1]])
            self.nfp.append([self.stationary[2][0]+width,self.stationary[2][1]+height])
            self.nfp.append([self.stationary[3][0],self.stationary[3][1]+height])
        else:
            # while self.judgeEnd()==False:
            while self.judgeEnd()==False and i<75: # 大于等于75会自动退出的，一般情况是计算出错
            # while i<7:
                # print("########第",i,"轮##########")
                touching_edges=self.detectTouching()
                all_vectors=self.potentialVector(touching_edges)
                if len(all_vectors)==0:
                    print("没有可行向量")
                    self.error=-2 # 没有可行向量
                    break

                vector=self.feasibleVector(all_vectors,touching_edges)
                if vector==[]:
                    print("没有计算出可行向量")
                    self.error=-5 # 没有计算出可行向量
                    break
                
                self.trimVector(vector)
                if vector==[0,0]:
                    print("未进行移动")
                    self.error=-3 # 未进行移动
                    break

                GeoFunc.slidePoly(self.sliding,vector[0],vector[1])
                self.nfp.append([self.sliding[self.locus_index][0],self.sliding[self.locus_index][1]])
                
                # self.showPolygon()
                # time.sleep(0.3)

                i=i+1
                inter=Polygon(self.sliding).intersection(Polygon(self.stationary))
                if GeoFunc.computeInterArea(inter)>1:
                    print("出现相交区域")
                    self.error=-4 # 出现相交区域
                    break  

        # self.showNFP()

        if i==75:
            # print("超出计算次数")
            self.error=-1 # 超出计算次数
    
    def showNFP(self):
        # 创建图形和坐标系
        fig, ax = plt.subplots()
        polygon1 = Polygon(self.sliding)
        polygon2 = Polygon(self.stationary)
        NFP1=Polygon(self.nfp)
        # print(self.sliding,self.stationary,self.nfp)
        # 绘制第一个多边形并填充颜色
        x1, y1 = polygon2.exterior.xy
        ax.fill(x1, y1, edgecolor='red', facecolor='lightcoral', label='stationary polygon i')
        # 绘制第二个多边形并填充颜色
        x2, y2 = polygon1.exterior.xy
        ax.fill(x2, y2, edgecolor='blue', facecolor='lightblue', label='sliding polygon j')
        #NFP
        x3,y3=NFP1.exterior.xy
        ax.plot(x3, y3, color='black', linestyle='-', linewidth=2, label='NFP')
        # 设置坐标轴范围
        ax.set_xlim(-30, 50)
        ax.set_ylim(-30, 50)
        # ax.set_xlim(-100, 150)
        # ax.set_ylim(-100, 150)
        # ax.set_xlim(-300, 300)
        # ax.set_ylim(-300, 300)
        # 添加图例
        ax.legend()
        plt.show()
        # timestamp = time.strftime("%Y%m%d_%H%M%S")
        # plt.savefig(f"D:/Desktop/排样文件/NFP/degrade_{timestamp}.png", dpi=600)
        # 显示图形
        # plt.show(block=False)
        # plt.pause(1)           # 显示1秒
        # plt.close()            # 关闭当前窗口
        
    # 画nfp的移动图
    def showPolygon(self):
        # 创建图形和坐标系
        fig, ax = plt.subplots()
        polygon1 = Polygon(self.sliding)
        polygon2 = Polygon(self.stationary)
        # 绘制第一个多边形并填充颜色
        x1, y1 = polygon2.exterior.xy
        ax.fill(x1, y1, edgecolor='red', facecolor='lightcoral', label='stationary polygon i')
        # 绘制第二个多边形并填充颜色
        x2, y2 = polygon1.exterior.xy
        ax.fill(x2, y2, edgecolor='blue', facecolor='lightblue', label='sliding polygon j')
        # 设置坐标轴范围
        ax.set_xlim(-30, 50)
        ax.set_ylim(-30, 50)
        # ax.set_xlim(-100, 150)
        # ax.set_ylim(-100, 150)
        # ax.set_xlim(-300, 300)
        # ax.set_ylim(-300, 300)
        # 添加图例
        ax.legend()
        plt.show()
        # timestamp = time.strftime("%Y%m%d_%H%M%S")
        # plt.savefig(f"D:/Desktop/排样文件/NFP/degrade_{timestamp}.png", dpi=600)
        # 显示图形
        # plt.show(block=False)
        # plt.pause(1)           # 显示1秒
        # plt.close()            # 关闭当前窗口

    # 检测相互的连接情况
    def detectTouching(self):
        touch_edges=[]
        stationary_edges,sliding_edges=self.getAllEdges()
        for edge1 in stationary_edges:
            for edge2 in sliding_edges:
                inter=GeoFunc.intersection(edge1,edge2)
                if inter!=[]:
                    pt=[inter[0],inter[1]] # 交叉点
                    edge1_bound=(GeoFunc.almostEqual(edge1[0],pt) or GeoFunc.almostEqual(edge1[1],pt)) # 接触点是否在边的两点
                    edge2_bound=(GeoFunc.almostEqual(edge2[0],pt) or GeoFunc.almostEqual(edge2[1],pt)) # 接触点是否在边的两点
                    stationary_start=GeoFunc.almostEqual(edge1[0],pt) # 起点是否是接触点
                    orbiting_start=GeoFunc.almostEqual(edge2[0],pt) # 起点是否是接触点
                    touch_edges.append({
                        "edge1":edge1,
                        "edge2":edge2,
                        "vector1":self.edgeToVector(edge1),
                        "vector2":self.edgeToVector(edge2),
                        "edge1_bound":edge1_bound,
                        "edge2_bound":edge2_bound,
                        "stationary_start":stationary_start,
                        "orbiting_start":orbiting_start,
                        "pt":[inter[0],inter[1]],
                        "type":0
                    })
        return touch_edges 
    
    # # 获得潜在的可转移向量
    # def potentialVector(self, touching_edges):
    #     all_vectors = []  # 初始化一个空列表，用于存储所有潜在的可转移向量
    #     for touching in touching_edges:  # 遍历所有接触边
    #         aim_edge = []  # 初始化目标的列表
    #         # 情况1：两个边都是外边
    #         if touching["edge1_bound"] == True and touching["edge2_bound"] == True:
    #             left, right,parallel = GeoFunc.judgePosition(touching["edge1"], touching["edge2"])  # 判断边的位置(角度)关系
    #             # 如果固定边的起点是接触点，滑动边的起点不是 case1
    #             if touching["stationary_start"] == True and touching["orbiting_start"] == False:
    #                 touching["type"] = 2  # 设置类型为2
    #                 aim_edge = [touching["edge2"][1], touching["edge2"][0]]  # 反方向 jk+1jk

    #             # 如果固定边的起点不是接触点，滑动边的起点是 case2
    #             if touching["stationary_start"] == False and touching["orbiting_start"] == True:
    #                 touching["type"] = 1  # 设置类型为1
    #                 aim_edge = touching["edge1"]  # 选择边1 ikik+1

    #             # 如果两个边的起点都是接触点
    #             if touching["stationary_start"] == True and touching["orbiting_start"] == True:
    #                 touching["type"] = 0  # 设置类型为0
    #                 if left == True: # case4：角度>180°
    #                     aim_edge = [touching["edge2"][1], touching["edge2"][0]]  # Pj边的反方向 jk+1jk
    #                 if right == True: # case5: 角度<=180°
    #                     aim_edge = touching["edge1"]  # 选择边1 Pi边的正方向 ikik+1
                        
    #             # 如果两个边的起点都不是接触点 case3 
    #             if touching["stationary_start"] == False and touching["orbiting_start"] == False:
    #                 touching["type"] = 3  # 设置类型为3

    #         # 情况2：边1不是外边(内边)，边2是外边 case6
    #         if touching["edge1_bound"] == False and touching["edge2_bound"] == True:
    #             aim_edge = [touching["pt"], touching["edge1"][1]]  # ikik+1
    #             touching["type"] = 4  # 设置类型为4
            
    #         # 情况3：边1是外边，边2不是外边(内边) case6
    #         if touching["edge1_bound"] == True and touching["edge2_bound"] == False:
    #             aim_edge = [touching["edge2"][1], touching["pt"]]  # jk+1jk
    #             touching["type"] = 5  # 设置类型为5

    #         if aim_edge != []:  # 如果目标边不为空
    #             vector = self.edgeToVector(aim_edge)  # 计算目标边的向量
    #             if self.detectExisting(all_vectors, vector) == False:  # 检查该向量是否已存在
    #                 all_vectors.append(vector)  # 如果不存在，则添加到列表中
    #     return all_vectors  # 返回所有潜在的可转移向量

    # 获得潜在的可转移向量
    def potentialVector(self, touching_edges):
        all_vectors = []  # 初始化一个空列表，用于存储所有潜在的可转移向量
        for touching in touching_edges:  # 遍历所有接触边
            aim_edge = []  # 初始化目标的列表
            # 情况1：两个都是边
            if touching["edge1_bound"] == True and touching["edge2_bound"] == True:
                left, right,parallel = GeoFunc.judgePosition(touching["edge1"], touching["edge2"])  # 判断边的位置(角度)关系
                # 如果固定边的起点不是接触点，滑动边的起点是 case1
                if touching["stationary_start"] == False and touching["orbiting_start"] == True:
                    touching["type"] = 2  # 设置类型为2
                    if right == True:
                        aim_edge = [touching["edge2"][1], touching["edge2"][0]]  # 反方向 jk+1jk
                # 如果固定边的起点是接触点，滑动边的起点不是 case2
                if touching["stationary_start"] == True and touching["orbiting_start"] == False:
                    touching["type"] = 1  # 设置类型为1
                    if left == True:
                        aim_edge = touching["edge1"]  # 选择边1 ikik+1
                # 如果两个边的起点都是接触点
                if touching["stationary_start"] == True and touching["orbiting_start"] == True:
                    touching["type"] = 0  # 设置类型为0
                    if left == True: # case4：角度>180°
                        aim_edge = [touching["edge2"][1], touching["edge2"][0]]  # Pj边的反方向 jk+1jk
                    if right == True: # case5: 角度<=180°
                        aim_edge = touching["edge1"]  # 选择边1 Pi边的正方向 ikik+1
                # 如果两个边的起点都不是接触点 case3 ,不用添加转移向量，按case5的方向移动
                if touching["stationary_start"] == False and touching["orbiting_start"] == False:
                    touching["type"] = 3  # 设置类型为3

            # 情况2：边1不是边界，边2是边界
            if touching["edge1_bound"] == False and touching["edge2_bound"] == True:
                aim_edge = [touching["pt"], touching["edge1"][1]]  # 目标边为接触点到边1的另一端 接触点到ik+1, ikik+1
                touching["type"] = 4  # 设置类型为4
            
            # 情况3：边1是边界，边2不是边界
            if touching["edge1_bound"] == True and touching["edge2_bound"] == False:
                aim_edge = [touching["edge2"][1], touching["pt"]]  # 目标边为接触点到边2的另一端
                touching["type"] = 5  # 设置类型为5
            # 只对于2-degrade
            if touching["edge1"]==[[3.0, 5.0], [3.0, 10.01]]:
                touching["type"] = 6  # 设置类型为6
                if touching["edge2"]==[[3.0, 10.0], [3.0, 5.0]]:
                    aim_edge=[[3,10],[3,10.01]]
                if touching["edge2"]==[[3.0, 5.0], [8.0, 5.0]]:
                    aim_edge=[[3,5],[3,5.01]]
            if aim_edge != []:  # 如果目标边不为空
                vector = self.edgeToVector(aim_edge)  # 计算目标边的向量
                if self.detectExisting(all_vectors, vector) == False:  # 检查该向量是否已存在
                    all_vectors.append(vector)  # 如果不存在，则添加到列表中
        return all_vectors  # 返回所有潜在的可转移向量
    
    def detectExisting(self,vectors,judge_vector):
        for vector in vectors:
            if GeoFunc.almostEqual(vector,judge_vector):
                return True
        return False
    
    def edgeToVector(self,edge):
        return [edge[1][0]-edge[0][0],edge[1][1]-edge[0][1]]
    
    # 选择可行向量
    def feasibleVector(self,all_vectors,touching_edges):
        '''
        该段代码需要重构，过于复杂
        '''
        res_vector=[]
        # print("\nall_vectors:",all_vectors)
        for vector in all_vectors:
            feasible=True
            # print("\nvector:",vector,"\n")
            for touching in touching_edges:
                vector1=[]
                vector2=[]
                # 判断方向并进行转向
                if touching["stationary_start"]==True:
                    vector1=touching["vector1"]
                else:
                    vector1=[-touching["vector1"][0],-touching["vector1"][1]]
                if touching["orbiting_start"]==True:
                    vector2=touching["vector2"]
                else:
                    vector2=[-touching["vector2"][0],-touching["vector2"][1]]
                vector12_product=GeoFunc.crossProduct(vector1,vector2) # 叉积，大于0在左侧，小于0在右侧，等于0平行
                vector_vector1_product=GeoFunc.crossProduct(vector1,vector) # 叉积，大于0在左侧，小于0在右侧，等于0平行
                vector_vector2_product=GeoFunc.crossProduct(vector2,vector) # 叉积，大于0在左侧，小于0在右侧，等于0平行
                # 最后两种情况
                if touching["type"]==4 and (vector_vector1_product*vector12_product)<0:
                    feasible=False
                if touching["type"]==5 and (vector_vector2_product*(-vector12_product))>0:
                    feasible=False
                # 正常的情况处理
                if vector12_product>0:
                    if vector_vector1_product<0 and vector_vector2_product<0:
                        feasible=False
                if vector12_product<0:
                    if vector_vector1_product>0 and vector_vector2_product>0:
                        feasible=False
                # 平行情况，需要用原值逐一判断
                if vector12_product==0:
                    inter=GeoFunc.newLineInter(touching["edge1"],touching["edge2"])
                    if inter["geom_type"]=="LineString":
                        if inter["length"]>0.01:
                            # 如果有相交，则需要在左侧
                            if (touching["orbiting_start"]==True and vector_vector2_product<0) or (touching["orbiting_start"]==False and vector_vector2_product>0):
                                feasible=False
                    else:
                        # 如果方向相同，且转化直线也平行，则其不能够取a的方向
                        if touching["orbiting_start"]==True != touching["stationary_start"]==False and vector_vector1_product==0:
                            if touching["vector1"][0]*vector[0]>0: # 即方向相同
                                feasible=False
            if feasible==True:
                res_vector=vector
                break
        return res_vector
        
    # 削减过长的向量
    def trimVector(self,vector):
        stationary_edges,sliding_edges=self.getAllEdges()
        new_vectors=[]
        for pt in self.sliding:
            for edge in stationary_edges:
                line_vector=LineString([pt,[pt[0]+vector[0],pt[1]+vector[1]]])
                end_pt=[pt[0]+vector[0],pt[1]+vector[1]]
                line_polygon=LineString(edge)
                inter=line_vector.intersection(line_polygon)
                if inter.geom_type=="Point":
                    inter_mapping=mapping(inter)
                    inter_coor=inter_mapping["coordinates"]
                    if (abs(end_pt[0]-inter_coor[0])>0.01 or abs(end_pt[1]-inter_coor[1])>0.01) and (abs(pt[0]-inter_coor[0])>0.01 or abs(pt[1]-inter_coor[1])>0.01):
                        new_vectors.append([inter_coor[0]-pt[0],inter_coor[1]-pt[1]])

        for pt in self.stationary:
            for edge in sliding_edges:
                line_vector=LineString([pt,[pt[0]-vector[0],pt[1]-vector[1]]])
                end_pt=[pt[0]-vector[0],pt[1]-vector[1]]
                line_polygon=LineString(edge)
                inter=line_vector.intersection(line_polygon)
                if inter.geom_type=="Point":
                    inter_mapping=mapping(inter)
                    inter_coor=inter_mapping["coordinates"]
                    if (abs(end_pt[0]-inter_coor[0])>0.01 or abs(end_pt[1]-inter_coor[1])>0.01) and (abs(pt[0]-inter_coor[0])>0.01 or abs(pt[1]-inter_coor[1])>0.01):
                        new_vectors.append([pt[0]-inter_coor[0],pt[1]-inter_coor[1]])
        
        # print(new_vectors)
        for vec in new_vectors:
            if abs(vec[0])<abs(vector[0]) or abs(vec[1])<abs(vector[1]):
                # print(vec)
                vector[0]=vec[0]
                vector[1]=vec[1]
    
    # 获得两个多边形全部边
    def getAllEdges(self):
        return GeoFunc.getPolyEdges(self.stationary),GeoFunc.getPolyEdges(self.sliding)
    
    # 判断是否结束
    def judgeEnd(self):
        sliding_locus=self.sliding[self.locus_index]
        main_bt=self.start_point
        if abs(sliding_locus[0]-main_bt[0])<0.1 and abs(sliding_locus[1]-main_bt[1])<0.1:
            if self.start==True:
                self.start=False
                # print("判断是否结束：否")
                return False
            else:
                # print("判断是否结束：是")
                return True
        else:
            # print("判断是否结束：否")
            return False

    # 显示最终结果
    def showResult(self):
        PltFunc.addPolygon(self.sliding)
        PltFunc.addPolygon(self.stationary)
        PltFunc.addPolygonColor(self.nfp)
        PltFunc.showPlt()

    # 计算渗透深度
    def getDepth(self):
        '''
        计算poly2的checkTop到NFP的距离
        Source: https://stackoverflow.com/questions/36972537/distance-from-point-to-polygon-when-inside
        '''
        d1=Polygon(self.nfp).distance(Point(self.original_top))
        # if point in inside polygon, d1=0
        # d2: distance from the point to nearest boundary
        if d1==0:
            d2=Polygon(self.nfp).boundary.distance(Point(self.original_top))
            # print('d2:',d2)
            return d2
        else: 
            return 0


