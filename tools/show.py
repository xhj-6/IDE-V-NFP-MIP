import matplotlib.pyplot as plt

class PltFunc(object):

    def addPolygon(poly):
        plt.fill(*zip(*poly), color="blue", alpha=0.5)
        for i in range(0,len(poly)):
            if i == len(poly)-1:
                PltFunc.addLine([poly[i],poly[0]])
            else:
                PltFunc.addLine([poly[i],poly[i+1]])

    def addPolygonColor(poly):
        
        for i in range(0,len(poly)):
            if i == len(poly)-1:
                PltFunc.addLine([poly[i],poly[0]],color="blue")
            else:
                PltFunc.addLine([poly[i],poly[i+1]],color="blue")

    def addLine(line,**kw):
        if len(kw)==0:
            plt.plot([line[0][0],line[1][0]],[line[0][1],line[1][1]],color="blue",linewidth=0.5)
        else:
            plt.plot([line[0][0],line[1][0]],[line[0][1],line[1][1]],color="blue",linewidth=0.5)            
    
    def showPlt(**kw):
        if len(kw)>0:
            if "minus" in kw:
                plt.axhline(y=0,c="blue")
                plt.axvline(x=0,c="blue")
                # plt.axis([-kw["minus"],kw["width"],-kw["minus"],kw["height"]])
                plt.axis([0,kw["width"],0,kw["length"]])
            else:
                plt.axis([0,kw["width"],0,kw["length"]])
        else:
            plt.axis([0,200,0,38])
            # plt.axis([-1000,2000,-979400.4498015114,20000])
            # plt.axis([-500,1000,0,1500])
        plt.show(block=False)
        plt.pause(1)           # 显示1秒
        plt.close()            # 关闭当前窗口
        # plt.show()
        # plt.clf()

    def showPolys(polys):
        for poly in polys:
            PltFunc.addPolygon(poly)
        PltFunc.showPlt(width=2000,height=2000)

    def saveFig(name):
        plt.savefig('figs\\'+name+'.png')
        plt.cla()
    
    