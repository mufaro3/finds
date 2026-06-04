from pprint import pprint
from .fish import *

def main():
    system = generate_system(
        n=3,
        bounds=[10,10,10],
        angle_delta=0.01
    )
    pprint(system)
    
if __name__ == '__main__':
    main()
