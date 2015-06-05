import sys
from ROOT import TLatex, TLegend


def myText(x,y,text, tsize,color):
    l = TLatex()
    l.SetTextSize(tsize); 
    l.SetNDC();
    l.SetTextColor(color);
    l.DrawLatex(x,y,text);


def getLegend(x1,y1,x2,y2,h1,h2,txt1,txt2,style1,style2):
    l = TLegend(x1,y1,x2,y2)
    l.SetFillColor(0)
    l.SetFillStyle(0)
    l.SetBorderSize(0)
    l.AddEntry(h1,txt1,style1)
    l.AddEntry(h2,txt2,style2)
    return l

    
