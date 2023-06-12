import re
class Dim2(object):
    def __init__(self, x_1, x_2):
        '''Defines x1 and x2 variables'''
        self.X_1 = x_1
        self.X_2 = x_2

    def getX1(self):
        return self.X_1

    def getX2(self):
        return self.X_2

    def __str__(self):
        return "(%s,%s)" % (self.X_1, self.X_2)



class Task(object):

    def __int__(self, start_point: Dim2, constraints: list, error: float, points: list):
        '''Define task with start conditions'''
        self.StartPoint = start_point
        self.Cond = constraints
        self.Error = error
        self.Points = points

    #getter
    @property
    def start_point(self):
        return self.StartPoint

    @property
    def conditions(self):
        return self.Cond

    @property
    def error(self):
        return self.Error

    @property
    def points(self):
        return self.Points


    '''Condition analyser'''










#def get_constraints():
    #n_con = int(input('Enter the number of constraints: '))
    #constraints = []
    #for i in range(0, n_con):
    #    constraint =





