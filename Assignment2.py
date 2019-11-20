import numpy as np
import functools
from matplotlib import pyplot as plt
####################################################
d = 3
print('type d is: ',type(d))
f = 4.544
print('type f is: ',type(f))
myName = 'Hassan'
print('type %s is: ' %myName, type(myName))
####################################################
myFavoriteFriuts = ['Apple','Orange','Peach','Banana']
for i in myFavoriteFriuts:
    print('My %d favorite fruit is' % myFavoriteFriuts.index(i),'\n')
####################################################
name1 = 'Hassan'
name2 = 'Ehsan'
brothers = '{} and {} are brothers.'.format(name1,name2)
print(brothers)
####################################################
myPhrase = 'I want {}'
print(myPhrase.format('an apple'))
####################################################
MySeris = map(lambda x: x**2, range(10))
print(list(MySeris),'\n')
####################################################
myListOfSchools = ['Harvard','MIT','Berkeley']
myListOfSchools.append('Stanford')
myListOfSchools.remove('MIT')
print('I like %d schools' %len(myListOfSchools))
####################################################
x = 1
y = 2
x,y = y,x
print('x,y =',x,',',y)
####################################################
myString1 = 'Hello'
myString2 = 'World.'
myString = myString1+' '+myString2
print('myString=', myString)
####################################################
dic = {}
dic['USA'] = 'Washington'
dic['Canada'] = 'Ottawa'
dic['Sweden'] = 'Stockholm'
for key in dic.keys():
    print('the capital of ' , key, 'is',dic[key])
####################################################
myPiNumber = 3.14
print(str(myPiNumber))
####################################################
import os
os.chdir(r'C:\Users\bahal\Desktop')
grades = open('Numbers','w')
grades.write('Hello,\n')
####################################################
#matrix = np.array([[1,2],[40,30]])
#matrix = np.zeros((4, 4))
#matrix = np.ones((4, 4))
matrix = np.random.uniform(0,1,(4, 4))
for (i, j), value in np.ndenumerate(matrix):
    print('value of: ',i+1,' and ',j+1,'=','%0.3f' %value)
####################################################
def Add(x_in,y_in):
    S = x_in + y_in
    return S
MyAdd = functools.reduce(Add,[1,2,3.55],2)
print('Sum of inputs is: %5.2f' %MyAdd)
####################################################
def GT_5 (x):
    if x>5: return True
    else: return False

MyFilter = filter(GT_5,range(11))
print('MyFilter is', list(MyFilter))
####################################################
class Student:
    numberOfStudents = 0
    def __init__(self, first, last, level, grade):
        self.first = first
        self.last = last
        self.level = level
        self.grade = grade
        Student.numberOfStudents += 1
    def fullname(self):
        return ('{}{}'.format(self.first, self.last))

    @classmethod
    def fromString(cls, studentString):
        first, last, level, grade = studentString.split('-')
        return cls(first, last, level, grade)

student1 = Student.fromString('Hassan-Bahaloo-12-A')
print('student name is: ', student1.first)
####################################################
class boneStrengthPlot:
    def __init__(self, rate, measure, abscissa, exponent):
        self.rate = rate
        self.measure = measure
        self.abscissa = abscissa
        self.exponent = exponent
    def makeRelation(self):
        strength = self.abscissa*self.measure**self.exponent*(rate/0.005)**0.06


####################################################
import pandas as pd
import numpy as np
np.random.seed(24)
df = pd.DataFrame({'A': np.linspace(1, 10, 10)})
df = pd.concat([df, pd.DataFrame(np.random.randn(10, 4), columns=list('BCDE'))],
               axis=1)
df.iloc[0, 2] = np.nan
print(df)
####################################################
count = 0
Nflip = 5000
theta = 0.5
HP = 0
head_Por = []
A = []
while count < Nflip:
    rand_Num = np.random.random()
    A.append(rand_Num)
    if rand_Num < theta:
       HP +=1
    count += 1
    head_Por.append(HP/count)

print('Length of head_Pro',len(head_Por))
plt.plot(range(0,count),head_Por)
plt.plot([0,Nflip],[theta,theta],'--')
plt.show()

#for i in range(0,Len(A))


