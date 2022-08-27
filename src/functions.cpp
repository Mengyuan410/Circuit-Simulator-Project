#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <complex>
#include <iterator>
#include <cmath>

double getnumber(std::string str)
{
    int i = 0;
    //extract numbers from string
    for (; i < str.length(); i++) { if (isdigit(str[i])) break; }

    str = str.substr(i, str.length() - i);
    //change string to double
    double number = atof(str.c_str());

    
    return number;
}

double get_value(std::string value) {
    
    std::string unit;
    double nvalue;
    unit = value[value.size() - 1];
    double number = getnumber(value);
   
    if (unit == "p") {
        nvalue = number * pow(10, -12);
    }
    else if (unit == "n") {
        nvalue = number * pow(10, -9);
    }
    else if (unit == "u") {
        nvalue = number * pow(10, -6);
    }
    else if (unit == "m") {
        nvalue = number * pow(10, -3);
    }
    else if (unit == "k") {
        nvalue = number * pow(10, 3);
    }
    else if (unit == "g") {
        nvalue = number * pow(10, 6);
    }
    else if (unit == "G") {
        nvalue = number * pow(10, 9);
    }
    else {
        nvalue = number;
    }
 
    return nvalue;
}

int get_nodes(std::string node) {
    
    std::string n = node;
    
    if (n == "0") {
        return 0;
    }
    else {
       
        std::string x = n.erase(0, 1);
        int intnode = std::stoi(x);
        
        return intnode;
    }
}

