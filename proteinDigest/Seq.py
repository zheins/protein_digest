class Seq:
    def __init__(self,seq,name):
        self.seq = seq
        self.name = name
    def length(self):
        return len(self.seq)
    def freq(self,sym):
        return self.seq.count(sym)
    def is_valid(self):
        for sym in self.seq:
            if sym not in valid_sym:
                return False
            return True
    
class Protein(Seq):
    mw = {'A':71.04,'C':103.01,'D':115.03,'E':129.04,'F':147.07,'G':57.02,'H':137.06,'I':113.08,'K':128.09,'L':113.08,'M':131.04,'N':114.04,'P':97.05,'Q':128.06,'R':156.10,'S':87.03,'T':101.05,'V':99.07,'W':186.08,'Y':163.06}
    valid_sym = 'ACDEFGHIKLMNPQRSTVWY'
    def molWt(self):
        return sum(map(self.mw.get,self.seq))
