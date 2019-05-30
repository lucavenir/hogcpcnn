from inout import *
from lib import *

def main():
    inputs = load_input('mysmallgraph')
    subgraph = [0, 1]
    print(u_triangle_count(inputs['adj']))

if __name__ == "__main__":
    main()
