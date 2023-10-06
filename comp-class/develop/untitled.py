def myFunction():
    print('this is a function')
    
myFunction()

def factorial(n):
    if n<0:
        print('undefined')
    else:
        f = 1
        for i in range(f,n+1):
            f = f*i
        return(f)
    
factorial(-4)