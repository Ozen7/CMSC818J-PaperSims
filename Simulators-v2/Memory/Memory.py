
class Memory:
    
    def __init__(self,bpc) -> None:
        self.bytespercycle = bpc
        self.remainding = bpc
    def request(self,numbytes):
        if 