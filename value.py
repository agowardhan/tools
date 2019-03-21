
class value: 
     
    '''
    Structure to save a value, error, fractional error, and unit. 
    '''
    
    def __init__(self, val=None, err=None, unit=None):
        
        self.value = val
        self.error = err
        if self.value and self.error: 
            self.frac  = self.error/self.value
        self.unit = unit 
    
    def calcFrac(self): 
        if self.value and self.error: 
            self.frac  = self.error/self.value
    
    def mult(self, const, unit=None): 
        
        self.value *= const
        self.error *= const
        if unit: 
            self.unit = unit 
            
            