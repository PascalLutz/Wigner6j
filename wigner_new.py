import math
import numpy as np
from pair_multiplication import *


class Fraction:
  """
  This class handles calculation without carrying out divisions or roots. It supports most basic operations and comes with a toString function.
  """
  def __init__(self, nominator=0, denominator=1, nominatorRoot = 1, denominatorRoot = 1):
    if isinstance(nominator, int):
      self.nominator = nominator
    else:
      self.nominator = round(nominator)
    if isinstance(denominator, int):
      self.denominator = denominator
    else:
      self.denominator = round(denominator)
    if isinstance(nominatorRoot, int):
      self.nominatorRoot = nominatorRoot
    else:
      self.nominatorRoot = round(nominatorRoot)
    if isinstance(denominatorRoot, int):
      self.denominatorRoot = denominatorRoot
    else:
      self.denominatorRoot = round(denominatorRoot)

  def __setattr__(self, name, value):
    if not isinstance(value, int):
      raise TypeError("value must be an int")
    super().__setattr__(name, value)

  def __nodenomroot__(self):
    return Fraction(self.nominator,  self.denominator * self.denominatorRoot, self.nominatorRoot * self.denominatorRoot, 1)

  def __add__(self, other):
    if isinstance(other, Fraction):
      if(self.denominatorRoot != 1):
        self = self.__nodenomroot__()
      if(other.denominatorRoot != 1):
        other = other.__nodenomroot__()
      if(self.nominatorRoot == other.nominatorRoot and self.denominatorRoot == 1 and other.denominatorRoot == 1):
        return Fraction(self.nominator * other.denominator + other.nominator * self.denominator, self.denominator * other.denominator, self.nominatorRoot,1).reduce()
      if(self.nominatorRoot == other.denominatorRoot and self.denominatorRoot == 1 and other.nominatorRoot == 1):
        return Fraction(self.nominator * other.denominator * other.denominatorRoot + other.nominator * self.denominator, self.denominator * other.denominator * other.denominatorRoot, self.nominatorRoot,1).reduce()
      if(self.nominatorRoot != 1 or self.denominatorRoot != 1 or other.nominatorRoot != 1 or other.denominatorRoot != 1):
        if(self==0):
          return other
        elif(other==0):
          return self
        if(self.nominatorRoot == 1 and other.nominatorRoot==1 and self.denominatorRoot == other.denominatorRoot):
          return Fraction(self.nominator*other.denominator + self.denominator*other.nominator,self.denominator*other.denominator,1,self.denominatorRoot).reduce()
        if(self.nominator == -other.nominator and self.denominator == other.denominator and self.nominatorRoot == other.nominatorRoot and self.denominatorRoot == other.denominatorRoot):
          return 0
        print("Adding two roots. Imprecision!", self, other)
        raise ValueError("Imprecision")
        return (self.nominator * other.denominator * math.sqrt(self.nominatorRoot * other.denominatorRoot) + other.nominator * self.denominator * math.sqrt(other.nominatorRoot * self.denominatorRoot))/(self.nominator*other.nominator*math.sqrt(self.denominatorRoot*other.denominatorRoot))
      return Fraction(self.nominator * other.denominator + other.nominator * self.denominator, self.denominator * other.denominator).reduce()
    elif isinstance(other, int):
      if(other==0):
        return self
      if(self.nominatorRoot != 1 or self.denominatorRoot != 1):
        print("Adding a int to a root. Imprecision!", self, other)
        return (other * self.denominator * math.sqrt(self.denominatorRoot) + self.nominator * math.sqrt(self.nominatorRoot)) / (self.denominator * math.sqrt(self.denominatorRoot))
      return Fraction(other * self.denominator + self.nominator, self.denominator).reduce()
    elif isinstance(other, float):
      return (other * self.denominator * math.sqrt(self.denominatorRoot) + self.nominator * math.sqrt(self.nominatorRoot)) / (self.denominator * math.sqrt(self.denominatorRoot))

  def __radd__(self, other):
    return self.__add__(other)

  def __sub__(self, other):# TODO with Roots
    if isinstance(other, Fraction):
      return Fraction(self.nominator * other.denominator - other.nominator * self.denominator, self.denominator * other.denominator).reduce()
    elif isinstance(other, int):
      return Fraction(other * self.denominator - self.nominator, self.denominator).reduce()
    elif isinstance(other, float):
      return Fraction(other * self.denominator - self.nominator, self.denominator).reduce()

  def __rsub__(self, other):
    return self.__sub__(other)

  def __mul__(self, other):
    if isinstance(other, Fraction):
      return Fraction(self.nominator * other.nominator * self.nominatorRoot * other.nominatorRoot,self.denominator * other.denominator, 1, self.nominatorRoot * self.denominatorRoot * other.nominatorRoot * other.denominatorRoot).reduce()
    elif isinstance(other, int):
      return Fraction(other * self.nominator, self.denominator, self.nominatorRoot, self.denominatorRoot).reduce()

  def __rmul__(self, other):
    return self.__mul__(other)

  #def __div__(self, other):
  #  print("BBBBBBBBBBB", self, other)
  #  if isinstance(other, Fraction):
  #    return Fraction(self.nominator * other.denominator * self.nominatorRoot * other.denominatorRoot,self.denominator * other.nominator, 1, self.nominatorRoot * self.denominatorRoot * other.nominatorRoot * other.denominatorRoot).reduce()
  #  elif isinstance(other, int):
  #    return Fraction(self.denominator*other, self.nominator, self.denominatorRoot, self.nominatorRoot).reduce()

  def __rtruediv__(self, other):
    if isinstance(other, Fraction):
      return Fraction(self.nominator * other.denominator * self.nominatorRoot * other.denominatorRoot,self.denominator * other.nominator, 1, self.nominatorRoot * self.denominatorRoot * other.nominatorRoot * other.denominatorRoot).reduce()
    elif isinstance(other, int):
      return Fraction(self.denominator*other, self.nominator, self.denominatorRoot, self.nominatorRoot).reduce()

  #def __future__(self, other):
  #  print("DDDDDDDDDDD", self, other)
  #  if(isinstance(self, int)):
  #    return Fraction(other.denominator*self, other.nominator, other.denominatorRoot, other.nominatorRoot)
  #  return self.__div__(other)
  
  def __truediv__(self, other):
    if isinstance(other, Fraction):
      return Fraction(self.nominator * other.denominator * self.nominatorRoot * other.denominatorRoot,self.denominator * other.nominator, 1, self.nominatorRoot * self.denominatorRoot * other.nominatorRoot * other.denominatorRoot).reduce()
    elif isinstance(other, int):
      return Fraction(self.nominator, other * self.denominator, self.nominatorRoot, self.denominatorRoot).reduce()

  def __lt__(self, other):
    if isinstance(other, Fraction):
      return (self.nominator * other.denominator * math.sqrt(self.nominatorRoot * other.denominatorRoot) < self.denominator * other.nominator * math.sqrt(self.denominatorRoot * other.nominatorRoot))
    elif isinstance(other, (int, float)):
      return self.nominator * math.sqrt(self.nominatorRoot) < self.denominator * other * math.sqrt(self.denominatorRoot)

  def __rlt__(self, other):
    return self.__lt__(other)

  def __eq__(self, other):
    if isinstance(other, Fraction):
      s = self.reduce()
      o = other.reduce()
      return (s.nominator == o.nominator) and (s.denominator == o.denominator) and (s.nominatorRoot == o.nominatorRoot) and (s.denominatorRoot == o.denominatorRoot)
    if isinstance(other, (float,int)):
      return self.denominator * other * self.denominatorRoot == self.nominator *math.sqrt(self.nominatorRoot)
    else:
      return False

  def __gt__(self, other):
    if isinstance(other, Fraction):
      return (self.nominator * other.denominator * math.sqrt(self.nominatorRoot * other.denominatorRoot) > self.denominator * other.nominator * math.sqrt(self.denominatorRoot * other.nominatorRoot))
    elif isinstance(other, (int, float)):
      return self.nominator * math.sqrt(self.nominatorRoot) > self.denominator * other * math.sqrt(self.denominatorRoot)

  def __rgt__(self, other):
    return self.__gt__(other)

  def sqrt(self):
    return Fraction(1,1,round(self.nominator * math.sqrt(self.nominatorRoot)), round(self.denominator * math.sqrt(self.denominatorRoot)))

  def sq(self):
    if(isinstance(self, int)):
      return Fraction(self*self,1)
    return Fraction(self.nominator**2 * self.nominatorRoot,self.denominator**2 * self.denominatorRoot)

  def __pow__(self,other):
    return Fraction(self.nominator**other, self.denominator**other, self.nominatorRoot**other, self.denominatorRoot**other)

  def abs(self):
    if self.nominator < 0:
      self.nominator = round(self.nominator * -1)
    return self

  def reduce(self):
    if self.denominator == 0:
      return None
    if self.denominator <= 0:
      self.denominator *= (-1)
      self.nominator *= (-1)
    if self.nominatorRoot != 1:
      self.nominator *= self.nominatorRoot
      self.denominatorRoot *= self.nominatorRoot
      self.nominatorRoot = 1
    i=2
    while(i**2<=self.denominatorRoot):
      if(self.denominatorRoot % i**2 == 0):
        self.denominatorRoot = round(self.denominatorRoot / (i**2))
        self.denominator *= i
        i-=1
      i+=1
    gcd = math.gcd(abs(self.nominator), abs(self.denominator))
    if gcd > 0:
      self.nominator = round(self.nominator / gcd)
      self.denominator = round(self.denominator / gcd)
    gcd = math.gcd(abs(self.nominator), abs(self.denominatorRoot))
    if(gcd>1):
      self.nominatorRoot = round(self.nominatorRoot * gcd)
      self.nominator = round(self.nominator / gcd)
      self.denominatorRoot = round(self.denominatorRoot / gcd)
    return self

  def __str__(self):
    
    if self.nominator == 0 and self.denominator != 0:
      return str(0)
    output = ""
    if(self.nominatorRoot == 1 and self.nominator == 1):
      output = "1"
    if(self.nominator != 1):
      output = str(self.nominator)
    if(self.nominator != 1 and self.nominatorRoot != 1):
      output += "*"
    if(self.nominatorRoot != 1 and self.denominatorRoot != 1):
      output += "Sqrt(" + str(self.nominatorRoot) + "/" + str(self.denominatorRoot) + ")"
    if(self.nominatorRoot != 1 and self.denominatorRoot == 1):
      output += "Sqrt(" + str(self.nominatorRoot) + ")"
    if((self.denominatorRoot == 1 and self.denominator == 1 ) or (self.nominatorRoot != 1 and self.denominatorRoot != 1 and self.denominator == 1)):
      return output
    else: output += "/"
    if(self.denominator != 1 and self.denominatorRoot != 1 and self.nominatorRoot == 1):
      output += "("
    if(self.denominator != 1):
      output += str(self.denominator)
    if(self.denominator != 1 and self.denominatorRoot != 1 and self.nominatorRoot == 1):
      output += "*"
    if(self.denominatorRoot != 1 and self.nominatorRoot == 1):
      output += "Sqrt(" + str(self.denominatorRoot) + ")"
    if(self.denominator != 1 and self.denominatorRoot != 1 and self.nominatorRoot == 1):
      output += ")"
    return output

  def __repr__(self):
    if self.nominator == 0 and self.denominator != 0:
      return f"0"
    output = str(self.nominator)
    if self.nominatorRoot != 1 : 
      output += "Sqrt[" + str(self.nominatorRoot) + "]"
    output += "/" + str(self.denominator)
    if self.denominatorRoot != 1 : 
      output += "Sqrt[" + str(self.denominatorRoot) + "]"
    return f"{output}"


class Vertex:
  """
  Handles Vertices and their young diagrams.
  """
  
  def __init__(self, young_diagram_1=[0], young_diagram_2=[0], young_diagram_3=[0], vertex_number = 1):
    if type(young_diagram_1) == YoungDiagram:
      self.young_diagram_1 = young_diagram_1
    else:
      self.young_diagram_1 = YoungDiagram(young_diagram_1)
    if type(young_diagram_2) == YoungDiagram:
      self.young_diagram_2 = young_diagram_2
    else:
      self.young_diagram_2 = YoungDiagram(young_diagram_2)
    if type(young_diagram_3) == YoungDiagram:
      self.young_diagram_3 = young_diagram_3
    else:
      self.young_diagram_3 = YoungDiagram(young_diagram_3)
    self.vertex_number = vertex_number

  def __str__(self):
    return f"{"First Representation: " + str(self.young_diagram_1) + " \nSecond Representation: " + str(self.young_diagram_2) + 
            " \nThird Representation: " + str(self.young_diagram_3) + " \nVertex Number: " + str(self.vertex_number)}"
    
  def __repr__(self):
    return f"{"First Representation: " + str(self.young_diagram_1) + " \nSecond Representation: " + str(self.young_diagram_2) + 
            " \nThird Representation: " + str(self.young_diagram_3) + " \nVertex Number: " + str(self.vertex_number)}"

  def is_well_defined(self):
    """
    Checks whether the vertex with the current information is well defined/ exists.
    """
    temp_young_diagram = YoungDiagram(self.young_diagram_3.partition, barred=True).evaluate_for_Nc(3)
    ds = self.young_diagram_1 * self.young_diagram_2
    ds = ds.evaluate_for_Nc(3)
    for i in range(ds.elements.size):
      if temp_young_diagram == ds.elements[i]:
        if self.vertex_number <= ds.multiplicities[i]:
          return True
    
    return False


class Wigner:
  def __init__(self, young_diagram_1=[0], young_diagram_2=[0], young_diagram_3=[0], young_diagram_4=[0], young_diagram_5=[0], young_diagram_6=[0], vertex_number_1 = 1, vertex_number_2 = 1, vertex_number_3 = 1, vertex_number_4 = 1):
    self.young_diagram_1 = YoungDiagram(young_diagram_1)
    self.young_diagram_2 = YoungDiagram(young_diagram_2)
    self.young_diagram_3 = YoungDiagram(young_diagram_3)
    self.young_diagram_4 = YoungDiagram(young_diagram_4)
    self.young_diagram_5 = YoungDiagram(young_diagram_5)
    self.young_diagram_6 = YoungDiagram(young_diagram_6)
    self.vertex_1 = Vertex(self.young_diagram_2,conjugate_diagram(self.young_diagram_3),conjugate_diagram(self.young_diagram_4),vertex_number_1)
    self.vertex_2 = Vertex(self.young_diagram_3,conjugate_diagram(self.young_diagram_1),self.young_diagram_5,vertex_number_2)
    self.vertex_3 = Vertex(self.young_diagram_1,conjugate_diagram(self.young_diagram_2),self.young_diagram_6,vertex_number_3)
    self.vertex_4 = Vertex(self.young_diagram_4,conjugate_diagram(self.young_diagram_5),conjugate_diagram(self.young_diagram_6),vertex_number_4)
    self.value = None
    self.cases = None

  def is_well_defined(self):
    """
    Checks whether the Wigner 6j-symbol with the current information is well defined/ exists.
    """
    return self.vertex_1.is_well_defined() and self.vertex_2.is_well_defined() and self.vertex_3.is_well_defined() and self.vertex_4.is_well_defined()

  def identify_case(self):
    if cases is None:
      cases = []
      if is_quark(self.vertex_2.young_diagram_1.partition) and is_quark(self.vertex_3.young_diagram_3.partition):
        cases.append(0)
      if is_quark(self.vertex_4.young_diagram_1.partition) and is_quark(self.vertex_3.young_diagram_3.partition) and is_gluon(self.vertex_2.young_diagram_3.partition):
        cases.append(1)
      if is_quark(self.vertex_2.young_diagram_1.partition) and is_gluon(self.vertex_3.young_diagram_3.partition):
        cases.append(2)
      if is_gluon(self.vertex_2.young_diagram_1.partition) and is_quark(self.vertex_3.young_diagram_3.partition):
        cases.append(2)
      if is_gluon(self.vertex_2.young_diagram_1.partition) and is_gluon(self.vertex_3.young_diagram_3.partition):
        cases.append(3)
      if is_gluon(self.vertex_4.young_diagram_1.partition) and is_gluon(self.vertex_3.young_diagram_3.partition) and is_gluon(self.vertex_2.young_diagram_3.partition):
        cases.append(4)
    return cases

  def get_value(self):
    if self.value is None:
      return 0 #TODO Implement
    return self.value




def conjugate_diagram(young_diagram, n=3):
  return YoungDiagram(young_diagram.partition, barred=True).evaluate_for_Nc(n)

def create_fractions_array(size):
  """
  Creates an Array if size is a tuple of two numbers or a List if size is a number containing Fractions-elements with nominator 0 and denominator 1.
  """
  array = []
  if type(size) is tuple:
    for i in range(size[0]):
      array.append([])
      for j in range(size[1]):
        array[i].append(Fraction())
  if isinstance(size, int):
    for i in range(size):
      array.append(Fraction())
  return array

def get_young_diagrams(numbers, sum, max_number=None):
  """
  Generates all possible Young diagrams with a certain amount of numbers and given total sum. max_number can be provided to limit the maximum possible digit.
  This is used to generate Young diagrams.
  """
  if sum == 0:
    return [np.zeros(numbers, dtype=int)]
  elif numbers == 1:
    return [[sum]]
  start = sum
  if max_number:
    start = min(sum, max_number)
  young_diagrams = False
  for i in range(start, math.ceil(sum / numbers)-1, -1):
    if type(young_diagrams) == bool:
      young_diagrams = np.insert(get_young_diagrams(numbers - 1, sum - i, i), 0, i, axis=1)
    else:
      young_diagrams = np.append(young_diagrams,np.insert(get_young_diagrams(numbers - 1, sum - i, i), 0, i, axis=1),axis=0)
  return young_diagrams

def create_young_tableaus(n, to_sum = 10, max_amount = 1000, includeSinglet = False):
  """
  Creates all Young diagrams with an amount of boxes up to n + to_sum and return them in order in a list.
  """
  diagrams = False
  sum = 1
  if(includeSinglet):
    sum = 0
  created = 0
  while sum <= to_sum and created <= max_amount:
    if type(diagrams) == bool:
      diagrams = get_young_diagrams(n, sum)
    else:
      diagrams = np.append(diagrams, get_young_diagrams(n, sum), axis=0)
      sum += 1
    created += 1
  return [diagram.tolist() for diagram in diagrams[1:]]

def get_number_of_boxes(alpha):
  """
  Returns the total amount of boxes in a diagram alpha.
  """
  count = 0
  for i in alpha:
    count += i
  return count

def get_factor(row, column):
  """
  Calculates the Factor of a specific cell in a diagram.
  """
  return column - row

def get_hooklength(diagram, row, column):
  """
  Calculates the Hooklength of a given diagram in a specific cell. The diagram needs to be given in a list of integers denoting the lengths of the rows.
  Returns the Hooklength as an Integer.
  """
  legLength = 0
  rowIndex = row
  while rowIndex < len(diagram):
    if diagram[rowIndex] >= column:
      legLength += 1
      rowIndex += 1
    else:
      break

  return diagram[row - 1] - column + 1 + legLength

def get_dimension(diagram, n):
  """
  Calculates the dimension of the irreducible represention of SU(N) of a given diagram via the factors-over-hooks formula.
  """
  if(type(diagram)==bool):
    return 0
  dimension = 1
  for row in range(1, len(diagram) + 1):
    for column in range(1, diagram[row - 1] + 1):
      # print("Field: ",row,column," Factor: ",get_factor(row, column)," HookLength: ",get_hooklength(diagram, row, column))
      # print((n+get_factor(row, column))/get_hooklength(diagram, row, column))
      dimension *= (n + get_factor(row, column)) / get_hooklength(diagram, row, column)
  return round(dimension)

def add_box_to_diagram(diagram, rowToAdd):
  """
  Add a box to a diagram which increments the number in that row by one. Makes a copy of the given diagram so it is not changed.

  Returns the diagram if it is admissable.
  Returns False is the diagram is not admissable.
  """
  newDiagram = diagram.copy()
  newDiagram[rowToAdd - 1] += 1
  if rowToAdd > 1:
    if newDiagram[rowToAdd - 1] > newDiagram[rowToAdd - 2]:
      return False
  if rowToAdd == len(newDiagram):
    for i in range(len(newDiagram)):
      newDiagram[i] -= 1
  if(is_viable_diagram(newDiagram)):
    return newDiagram
  else: return False

def remove_box_from_diagram(diagram, rowToRemove):
  """
  Removes a box from a diagram. 

  Returns the diagram if it is admissable and False otherwise.
  """
  newDiagram = diagram.copy()
  if(newDiagram[rowToRemove - 1] == 0):
    for i in range(len(diagram)-1):
      newDiagram[i] += 1
    return newDiagram
  newDiagram[rowToRemove - 1] -= 1
  if rowToRemove > 1:
    if(rowToRemove<len(diagram)):
      if newDiagram[rowToRemove - 1] < newDiagram[rowToRemove]:
        return False
  if(is_viable_diagram(newDiagram)):
    return newDiagram
  else: return False

def can_add_box(diagram, rowToAdd):
  """
  Returns whether a box can be added to the diagram while still being admissable.
  """
  if rowToAdd > 1:
    if diagram[rowToAdd - 1] >= diagram[rowToAdd - 2]:
      return False
  return True

def calculate_dimensions(diagram, n):
  """
  Returns two variables containing all dimensions of the diagrams obtained by adding one or two boxes respectively.
  These are sorted by which row the boxes are added in.
  If any diagram is not admissable the dimension will be 0.
  """

  dimensionsOneBox = np.zeros(n, dtype=int)
  dimensionsTwoBoxes = np.zeros([n, n], dtype=int)

  for i in range(1, n + 1):
    if not can_add_box(diagram, i):
      continue
    else:
      dimensionsOneBox[i - 1] = get_dimension(add_box_to_diagram(diagram, i), n)

  for row1 in range(1, n + 1):
    for row2 in range(1, n + 1):
      if can_add_box(diagram, row1):
        if can_add_box(add_box_to_diagram(diagram, row1), row2):
          dimensionsTwoBoxes[row1 - 1][row2 - 1] = get_dimension(add_box_to_diagram(add_box_to_diagram(diagram, row1), row2), n)
      elif can_add_box(diagram, row2):
        if can_add_box(add_box_to_diagram(diagram, row2), row1):
          dimensionsTwoBoxes[row1 - 1][row2 - 1] = get_dimension(add_box_to_diagram(add_box_to_diagram(diagram, row2), row1), n)
      else:
        dimensionsTwoBoxes[row1 - 1][row2 - 1] = 0
  return dimensionsOneBox, dimensionsTwoBoxes

def calculate_A_matrix(diagram, n, dimensionsOneBox, dimensionsTwoBoxes):
  """
  Calculates the lefthandside matrix found in appendix B2 of the quark-paper.
  """
  A = create_fractions_array((n, n))
  for row1 in range(1, n + 1):
    for row2 in range(row1, n + 1):
      if row1 == row2:
        continue
      else:
        k = get_hooklength(diagram, row1, diagram[row2 - 1] + 1)
        if dimensionsOneBox[row1 - 1] == 0:
          value1 = Fraction()
        else:
          value1 = Fraction(dimensionsTwoBoxes[row1 - 1][row2 - 1],dimensionsOneBox[row1 - 1]) * Fraction(1, k + 1)
        A[row1 - 1][row2 - 1] = value1
        if dimensionsOneBox[row2 - 1] == 0:
          value2 = Fraction()
        else:
          value2 = (-1* value1* Fraction(dimensionsOneBox[row1 - 1], dimensionsOneBox[row2 - 1]))
        A[row2 - 1][row1 - 1] = value2
  return A

def calculate_rhs(dimensionsOneBox, dimensionsTwoBoxes, n):
  """
  Calculates the righthandside vector found in appendix B2 of the quark-paper.
  """
  rhs = create_fractions_array(n)
  for i in range(1, n + 1):
    if dimensionsOneBox[i - 1] == 0:
      value = Fraction()
    else:
      value = 1 - Fraction(dimensionsTwoBoxes[i - 1][i - 1], dimensionsOneBox[i - 1])
    rhs[i - 1] = value
  return rhs

def find_changed_index(alpha, beta):
  """
  Returns the index where alpha and beta are different. Returns n when a box was added in the last row.
  Returns -1 if beta can not be constructed by adding a box to alpha.
  """
  index = -1
  reduction = True
  equalIndex = -1
  firstOccurence = True
  singleBox = True
  n = len(alpha)
  counter = 0
  for i in range(n):
    if(abs(alpha[i]-beta[i])==1):
      counter += 1
      index = i
  if(counter==1):
    return index
  elif(counter==n-1):
    return n-1
  else: 
    return -1
    if(alpha[i] != beta[i]):
      if(alpha[i] == beta[i]+1):
        if(firstOccurence):
          firstOccurence = False
          index = i
        else: singleBox = False
      else: singleBox = False
  if(singleBox):
    return index
  if(reduction==True and alpha[n-1] == beta[n-1]):
    return n-1
  return -1

def find_intermediate_diagram(alpha, beta, index=1):
  """
  Returns the diagram from which both alpha and beta can be constructed via removing a box. 
  If alpha and beta are equal, index can be used to determine where a box shall be added to alpha.
  """
  diagrams = []
  diagrams1 = []
  diagrams2 = []
  n = len(alpha)
  for i in range(n):
    diagrams1.append(add_box_to_diagram(alpha,i+1))
    diagrams2.append(add_box_to_diagram(beta,i+1))
  for diagram in diagrams1:
    if diagram in diagrams2 and diagram not in diagrams:
      diagrams.append(diagram)
  #print(diagrams, alpha, beta, diagrams1, diagrams2)
  if(len(diagrams)>=index):
    return diagrams[index-1]
  else:
    return False
  
def find_intermediate_diagrams(alpha, beta, checkForChargeConjugation=True, index=1):
  """
  Returns all diagrams from which both alpha and beta can be constructed via removing a box. 
  If alpha and beta are equal, index can be used to determine where a box shall be added to alpha.
  """
  diagrams = []
  diagrams1 = []
  diagrams2 = []
  chargeConjugate = False
  n = len(alpha)
  #if(checkForChargeConjugation):
  #  if(get_number_of_boxes(alpha) >= get_number_of_boxes(conjugate(alpha, n)) and get_number_of_boxes(beta) >= get_number_of_boxes(conjugate(beta, n))):
  #    print("charge", alpha, beta)
  #    alpha = conjugate(alpha, n)
  #    beta = conjugate(beta, n)
  #    chargeConjugate = True
    #if(get_number_of_boxes(beta) > get_number_of_boxes(conjugate(beta, n))):
      
     # chargeConjugate = True
  for i in range(n):
    diagram = add_box_to_diagram(alpha,i+1)
    if(not type(diagram) == bool):
      diagrams1.append(diagram)
    diagram = remove_box_from_diagram(beta,i+1)
    #if(not type(diagram) == bool):
    #  diagrams2.append(diagram)
  for diagram in diagrams1:
    if((diagram not in diagrams) and (diagram != False)):
      for i in range(n):
        if(remove_box_from_diagram(diagram, i+1) == beta):
          if((diagram not in diagrams) and (diagram != False)):
            diagrams.append(diagram)
    #if ((diagram in diagrams2) and (diagram not in diagrams) and (diagram != False)):
    #  diagrams.append(diagram)
  if(checkForChargeConjugation):
    for i in range(0, len(diagrams)):
      if(chargeConjugate):
        diagrams[i] = conjugate(diagrams[i], n)
        print(diagrams)
  return diagrams

def find_child_diagram(alpha, beta):
  """
  Returns a list of all diagrams that are one box smaller than alpha and beta and can be constructed by removing a box from them. 
  The list contains only a single diagram if alpha and beta are different.
  Return False if no diagram is admissable.
  """
  diagrams = []
  diagrams1 = []
  diagrams2 = []
  n = len(alpha)
  for i in range(n):
    diagrams1.append(remove_box_from_diagram(alpha,i+1))
    diagrams2.append(remove_box_from_diagram(beta,i+1))
  for diagram in diagrams1:
    if diagram in diagrams2 and diagram not in diagrams and diagram != False:
      diagrams.append(diagram)
  diagrams.sort(reverse=True)
  return diagrams

def conjugate(alpha, n):
  """
  Returns the conjuagte diagrams
  """
  newDiagram = alpha.copy()
  for i in range(1,n-1):
    newDiagram[i] = newDiagram[0] - newDiagram[i]
  return newDiagram

def is_viable_diagram(alpha):
  """
  Checks wheater a diagram is a valid Young diagram.
  """
  for i in range(len(alpha)-1):
    if(alpha[i]<alpha[i+1]):
      return False
  return True

def is_quark(alpha):
  """
  Checks whether a diagram is a quark.
  """
  return (alpha == None or alpha == "3" or alpha == [] or alpha == [1,0,0])

def is_adjoint(alpha,n=3):
  """
  Checks whether a diagram is a the adjoint representation.
  """
  test = [n-1,1]
  for i in range(2,n):
    test.append(0)
  return alpha == test

def is_gluon(alpha, n=3):
  """
  Checks whether a diagram is a gluon.
  """
  return (alpha == "gluon" or is_adjoint(alpha, n))

def is_real(alpha): # TODO compatible with any n
  return alpha[0] == alpha[1] *2 and alpha[2] == 0

def kronecker_delta(i,j):
  """
  Returns 1 if both inputs are equal and 0 otherwise.
  """
  if i==j : return 1
  return 0

def scalar_product(j,k,alpha,n=3):
  """
  Returns the scalarproduct as defined in the gluon-paper.
  """
  product = Fraction(1,(n**2)-1)
  if(kronecker_delta(j,k) and find_intermediate_diagrams(alpha, alpha)[j-1]):
    #print("j",j,"k",k,Fraction(1,get_dimension(find_intermediate_diagrams(alpha, alpha)[j-1],n)),Fraction(1,n*get_dimension(alpha,n)),Fraction(1,get_dimension(find_intermediate_diagrams(alpha, alpha)[j-2],n)) - Fraction(1,n*get_dimension(alpha,n)))
    product *= Fraction(1,get_dimension(find_intermediate_diagrams(alpha, alpha)[j-1],n)) - Fraction(1,n*get_dimension(alpha,n))
  else : product *= Fraction(-1,n*get_dimension(alpha,n))
  return product

def calculate_coefficient(alpha, gamma, a, i, n, debug=False):
  """
  Returns the coefficients as defined in the gluon-paper.
  """
  #if(alpha == [0,0,0]): #TODO Compatible with n
  #  return 0
  if(a==1 and i == 2):
    return 0
  if(alpha != gamma):
    if(a == 2 or i == 2):
      return 0
    lambdas = find_intermediate_diagrams(alpha, gamma)
    if(len(lambdas) == 0 or type(lambdas) == bool):
      return 0
    lambda1 = lambdas[0]
    if(debug):
      print("Coefficient Calculation with lambda, alpha, gamma:",lambda1, alpha, gamma)
    if(lambda1 == False):
      return 0
    if(debug):
      #TODO Check if there even is a second vertex
      print("Coefficient:" , get_dimension(lambda1,n),Fraction(1,1,get_dimension(lambda1,n) * (n**2 -1),1).reduce())
    return Fraction(1,1,get_dimension(lambda1,n) * (n**2 -1),1).reduce()
  else:
    diagrams = find_intermediate_diagrams(alpha, alpha)
    if(len(diagrams)>1):
      if(a == len(diagrams)):
        return 0
    if(a == 1):
      return 1/Fraction.sqrt(scalar_product(1,1,alpha,n)).reduce()
    elif(a == 2):
      if(i == 1):
        if((scalar_product(1,1,alpha,n)*scalar_product(2,2,alpha,n)-scalar_product(1,2,alpha,n)**2) == 0):
          print("coefficient not existing")
          return 0
        return Fraction.sqrt(scalar_product(1,1,alpha,n)/(scalar_product(1,1,alpha,n)*scalar_product(2,2,alpha,n)-scalar_product(1,2,alpha,n)**2))*(-1*scalar_product(1,2,alpha,n)/scalar_product(1,1,alpha,n)).reduce()
      elif(i == 2):
        if((scalar_product(1,1,alpha,n)*scalar_product(2,2,alpha,n)-scalar_product(1,2,alpha,n)**2) == 0):
          print("coefficient not existing")
          return 0
        return Fraction.sqrt(scalar_product(1,1,alpha,n)/(scalar_product(1,1,alpha,n)*scalar_product(2,2,alpha,n)-scalar_product(1,2,alpha,n)**2)).reduce() 

def calculate_normalization_plus(alpha, alpha_prime, n, debug):
  norm = 1 / Fraction.sqrt(Fraction(2,(n**2-1)) * (Fraction(1,get_dimension(alpha_prime,n))+calculate_wigner_6j_two_quarks(alpha,alpha_prime, [1,0,0], alpha, conjugate(alpha_prime,n), [1,0,0], n, debug) - Fraction(2,n*get_dimension(alpha,n))))
  if debug:
    print(debug * "  "+"Normalization plus: ",norm, "= root of: ", Fraction(2,(n**2-1)) ,"*(", Fraction(1,get_dimension(alpha_prime,n)),"+",calculate_wigner_6j_two_quarks(alpha,alpha_prime, [1,0,0], alpha, conjugate(alpha_prime,n), [1,0,0], n, debug),"-",Fraction(4,n*get_dimension(alpha,n)),")")
  return norm

def calculate_normalization_minus(alpha, alpha_prime, n, debug):
  norm = 1 / Fraction.sqrt(Fraction(2,(n**2-1)) * (Fraction(1,get_dimension(alpha_prime,n))-calculate_wigner_6j_two_quarks(alpha,alpha_prime, [1,0,0], alpha, conjugate(alpha_prime,n), [1,0,0], n, debug)))
  if debug:
    print(debug * "  "+"Normalization minus: ",norm, "= root of: ", Fraction(2,(n**2-1)) , Fraction(1,get_dimension(alpha_prime,n)),calculate_wigner_6j_two_quarks(conjugate(alpha_prime,n),alpha,[1,0,0],alpha_prime,alpha,[1,0,0], n, debug))
  return norm

def check_vertex_sign(vertex, alpha, beta, gamma, delta, epsilon, zeta, v1,v2,v3,v4,n=3):
  if type(vertex) is int:
    match vertex:
      case 1:
        reps = [beta, gamma, delta,v1]
      case 2:
        reps = [alpha, gamma, epsilon,v2]
      case 3:
        reps = [alpha, beta, zeta,v3]
      case 4:
        reps = [delta, epsilon, zeta,v4]
    if ((reps[0] == [1,0,0] and reps[1] == [1,0,0] and reps[2] == [1,1,0] and reps[3]==1)
      or (reps[0] == [1,0,0] and reps[1] == [1,1,0] and reps[2] == [1,0,0] and reps[3]==1)
       or (reps[0] == [1,1,0] and reps[1] == [1,0,0] and reps[2] == [1,0,0] and reps[3]==1)):
       return -1
    elif ((reps[0] == [2,1,0] and reps[1] == [4,2,0] and reps[2] == [4,2,0] and reps[3]==1)
      or (reps[0] == [4,2,0] and reps[1] == [4,2,0] and reps[2] == [2,1,0] and reps[3]==1)
       or (reps[0] == [4,2,0] and reps[1] == [2,1,0] and reps[2] == [4,2,0] and reps[3]==1)):
       return -1
    else: return 1
  elif type(vertex) is list:
    sign = 1
    for value in vertex:
      sign *= check_vertex_sign(value, alpha, beta, gamma, delta, epsilon, zeta,v1,v2,v3,v4,n)
    return sign
  else: return 1

def young_addition(alpha, beta):
  """
  Adds the diagram obtained by adding rowlengths of both provided diagrams.
  """
  aCopy = alpha.copy()
  for i in range(len(alpha)):
    if(beta[i]):
      aCopy[i] += beta[i]
  return aCopy

def young_multiplication_quark(alpha):
  """
  Carries out young multiplication of alpha and a quark.
  """
  diagrams = []
  for i in range(1,len(alpha)+1):
    diag = add_box_to_diagram(alpha,i)
    if diag != False:
      diagrams.append(diag)
  return diagrams

def young_multiplication_antiquark(alpha):
  """
  Carries out young multiplication of alpha and a quark.
  """
  diagrams = [[1,1,0],[1,0,1],[0,1,1]]
  new_diagrams = []
  for diagram in diagrams:
    new_diagram = young_addition(diagram,alpha)
    if new_diagram != False:
      if new_diagram[len(new_diagram)-1] == 1:
        for i in range(len(new_diagram)):
          new_diagram[i] -= 1
      new_diagrams.append(new_diagram)
  return new_diagrams

def young_multiplication_gluon(alpha, unique=True):
  """
  Carries out young multiplication of alpha and a gluon.
  """
  diagrams = []
  intermediateDiagrams = []
  for i in range(1,len(alpha)+1):
    diag = add_box_to_diagram(alpha,i)
    if(diag != False):
      intermediateDiagrams.append(diag)
  for i in intermediateDiagrams:
    diag2 = find_child_diagram(i,i)
    for j in diag2:
      if j != False:
        if unique:
          if j not in diagrams:
            diagrams.append(j)
        else:
          diagrams.append(j)
  return diagrams

def young_multiplication(alpha, beta):
  """
  Carries out young multiplication with both diagrams.
  """
  #TODO compatible with n
  aCopy = alpha.copy()
  bCopy = beta.copy()
  n = len(aCopy)
  diagrams = [aCopy]
  for i in range(n):
    firstOrderDiagrams = []
    for a in diagrams:
      if(beta[i]==0):
        break
      additions = create_young_tableaus(n,beta[i])
      print(additions)
      for j in additions:
        possibleDiagram = young_addition(a, j)
        print("pD: ",possibleDiagram, "alpha ", aCopy)
        for k in range(1,n):
          if(possibleDiagram[k]>aCopy[k-1]):
            print("violation",k)
            continue
        firstOrderDiagrams.append(possibleDiagram)
    diagrams += firstOrderDiagrams
  consistencySum = 0
  for i in diagrams:
    consistencySum += get_dimension(i,3)
  print(consistencySum, get_dimension(aCopy,3)*get_dimension(bCopy,3))
  return diagrams

def calculateSiiii(alpha, beta, delta, epsilon, row1, n):
  """
  Calculates a specific 2 quark symbol.
  """
  return Fraction(1,get_dimension(alpha, n))

def calculateSijii(alpha, beta, delta, epsilon, row1, row2, n):
  """
  Calculates a specific 2 quark symbol.
  """
  return Fraction(-1,get_dimension(alpha, n)) * Fraction(1,get_hooklength(epsilon, row1, epsilon[row2 - 1] + 1)+1)

def calculateSijjj(alpha, beta, delta, epsilon, row1, row2, n):
  """
  Calculates a specific 2 quark symbol.
  """
  return Fraction(1,get_dimension(alpha, n)) * Fraction(1,get_hooklength(epsilon, row1, epsilon[row2 - 1] + 1)+1)

def calculateSijij(alpha, beta, delta, epsilon, row1, row2, n):
  """
  Calculates a specific 2 quark symbol.
  """
  return Fraction(1,1,1,get_dimension(epsilon,n)*get_dimension(beta,n)).reduce()

def calculate_wigner_6j_two_quarks(alpha, beta, gamma, delta, epsilon, zeta, v1=1, v2=1, v3=1, v4=1, n=3, debug=True):
  """
  Calculates a case 0 6j-symbol.
  """
  if(debug):
    print(debug * "  "+"TwoQuarkSymbol", alpha, beta, gamma, delta, epsilon, zeta)
  j = find_changed_index(epsilon, delta)
  i = find_changed_index(epsilon, alpha)
  j2 = find_changed_index(alpha, beta)
  i2 = find_changed_index(delta, beta)
  #if(debug):
    #print(i,j,i2,j2,alpha, beta, delta,epsilon)
  if(j==-1 or i==-1 or j2==-1 or i2==-1):
    return 0
  if(j==i and j==j2 and j==i2):
    return calculateSiiii(alpha, beta, delta, epsilon, j+1, n)
  if(j==i and i2==j2 and j<j2):
    return calculateSijii(alpha, beta, delta, epsilon, j+1, j2+1, n)
  if(j==i and i2==j2 and j>j2):
    return calculateSijjj(alpha, beta, delta, epsilon, j2+1, j+1, n)
  if(j==j2 and i==i2 and i!=j):
    if(i<j):
      return calculateSijij(alpha, beta, delta, epsilon, i+1, j+1, n)
    else:
      return calculateSijij(alpha, beta, delta, epsilon, j+1, i+1, n)
  return False

def is_6j_with_two_quarks(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug):
  """
  Checks whether the provided 6j is a case 0 symbol.
  """
  return(is_quark(gamma) and is_quark(zeta))

def is_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Checks whether the provided 6j is a case 1 symbol.
  """
  #inner vertex cases
  if(find_changed_index(alpha, beta) != -1 and is_gluon(epsilon) and is_quark(delta) and is_quark(zeta) and v1 == 1 and v3 == 1 and v4 == 1):
    return True # calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  #if(find_changed_index(beta, alpha) != -1 and is_gluon(delta) and is_quark(epsilon) and is_quark(conjugate(zeta, n)) and v2 == 1 and v3 == 1 and v4 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(gamma, alpha, beta, conjugate(zeta, n), delta, epsilon, v3, v1, v2, v4, n, debug)
  #TODO there are more cases
  if(is_gluon(epsilon) and is_quark(conjugate(delta, n)) and is_quark(conjugate(zeta, n)) and v1 == 1 and v2 == 1 and v4 == 1):
    return True # calculate_6j_with_quark_gluon_vertex(gamma, beta, alpha, [1,0,0], epsilon, [1,0,0], v3, v2, v1, v4, n, debug)
  #bottom right vertex cases]
  #if(find_changed_index(alpha, epsilon) != -1 and is_gluon(delta) and is_quark(beta) and is_quark(gamma) and v1 == 1 and v2 == 1 and v3 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(epsilon, alpha, conjugate(zeta, n), beta, delta, gamma, v3, v4, v2, v1, n, debug)
  #if(find_changed_index(alpha, zeta) != -1 and is_gluon(gamma) and is_quark(beta) and is_quark(delta) and v1 == 1 and v3 == 1 and v4 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(conjugate(alpha, n), zeta, conjugate(epsilon, n), delta, gamma, beta, v4, v2, v3, v1, n, debug)
  #if(find_changed_index(alpha, epsilon) != -1 and is_gluon(beta) and is_quark(conjugate(gamma, n)) and is_quark(delta) and v1 == 1 and v2 == 1 and v4 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(zeta, epsilon, alpha, conjugate(gamma, n), beta, delta, v2, v3, v4, v1, n, debug)
  #bottom left vertex cases
  #if(find_changed_index(epsilon, delta) != -1 and is_gluon(alpha) and is_quark(beta) and is_quark(zeta) and v1 == 1 and v3 == 1 and v4 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(epsilon, delta, conjugate(gamma, n), beta, alpha, zeta, v1, v2, v4, v3, n, debug)
  #if(find_changed_index(conjugate(epsilon, n), gamma) != -1 and is_gluon(zeta) and is_quark(beta) and is_quark(alpha) and v1 == 1 and v2 == 1 and v3 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(conjugate(delta, n), gamma, conjugate(epsilon, n), alpha, zeta, beta, v2, v4, v1, v3, n, debug)
  #if(find_changed_index(delta, epsilon) != -1 and is_gluon(beta) and is_quark(conjugate(zeta, n)) and is_quark(alpha) and v2 == 1 and v3 == 1 and v4 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(conjugate(gamma, n), epsilon, delta, conjugate(zeta, n), beta, alpha, v4, v1, v2, v3, n, debug)
  #top vertex cases
  #if(find_changed_index(delta, beta) != -1 and is_gluon(epsilon) and is_quark(alpha) and is_quark(gamma) and v1 == 1 and v2 == 1 and v3 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(zeta, beta, delta, gamma, epsilon, alpha, v1, v4, v3, v2, n, debug)
  #if(find_changed_index(beta, delta) != -1 and is_gluon(alpha) and is_quark(epsilon) and is_quark(conjugate(gamma, n)) and v1 == 1 and v2 == 1 and v4 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(beta, delta, zeta, epsilon, alpha, conjugate(gamma, n), v4, v3, v1, v2, n, debug)
  #if(find_changed_index(beta, zeta) != -1 and is_gluon(gamma) and is_quark(conjugate(alpha, n)) and is_quark(conjugate(epsilon, n)) and v2 == 1 and v3 == 1 and v4 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(delta, zeta, beta, conjugate(alpha, n), gamma, conjugate(epsilon, n), v3, v1, v4, v2, n, debug)
  return False

def is_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Checks whether the provided 6j is a case 2 symbol.
  """
  #Cases with Gluon as inner line
  #case 1:
  #if(is_quark(alpha) and is_gluon(delta) and v2 == 1 and v3 == 1 and (v4 == 1 or v4 == 2)):
  #  return calculate_6j_with_quark_gluon_opposing(beta, gamma, alpha, conjugate(epsilon, n), zeta, delta, v2, v3, v1, v4, n, debug)
  #case 2:
  #if(is_quark(beta) and is_gluon(epsilon) and v1 == 1 and v3 == 1 and (v4 == 1 or v4 == 2)):
  #  return calculate_6j_with_quark_gluon_opposing(gamma, alpha, beta, conjugate(zeta, n), conjugate(delta, n), epsilon, v3, v1, v2, v4, n, debug)
  #case 3:
  if(is_quark(gamma) and is_gluon(zeta) and v1 == 1 and v2 == 1):
    return True
    return calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  #Cases with Gluon as outer line
  #case 1:
  #if(is_gluon(alpha) and is_quark(delta) and v3 == 1 and v1 == 1 and (v2 == 1 or v2 == 2)):
  #  return calculate_6j_with_quark_gluon_opposing(zeta, beta, delta, gamma, conjugate(epsilon, n), alpha, v2, v3, v1, v4, n, debug)
  #case 2:
  #if(is_gluon(beta) and is_quark(epsilon) and v2 == 1 and v4 == 1 and (v1 == 1 or v1 == 2)):
  #  return calculate_6j_with_quark_gluon_opposing(alpha, conjugate(zeta, n), epsilon, conjugate(delta, n), gamma, beta, v4, v2, v3, v1, n, debug)
  #case 3:
  #if(is_gluon(gamma) and is_quark(zeta) and v2 == 1 and v4 == 1 and (v2 == 1 or v2 == 2)):
  #  return calculate_6j_with_quark_gluon_opposing(beta, delta, zeta, epsilon, alpha, gamma, v4, v3, v1, v2, n, debug)
  return False

def is_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Checks whether the provided 6j is a case 3 symbol.
  """
  #case 1:
  if(is_gluon(gamma) and is_gluon(zeta) and (v4 == 1 or v4 == 2)):
    return calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  #case 2:
  #if(is_gluon(beta) and is_gluon(epsilon) and (v4 == 1 or v4 == 2)):
  #  return calculate_6j_with_two_gluon(gamma, alpha, beta, conjugate(zeta, n), conjugate(delta, n), epsilon, v3, v1, v2, v4, n, debug)
  #case 3:
  #if(is_gluon(alpha) and is_gluon(delta) and (v4 == 1 or v4 == 2)):
  #  return calculate_6j_with_two_gluon(beta, gamma, alpha, conjugate(epsilon, n), zeta, alpha, v2, v3, v1, v4, n, debug)
  return False

def is_6j_with_three_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Checks whether the provided 6j is a case 4 symbol.
  """
  #correlates to vertex 4 being triple gluon
  if(is_gluon(delta) and is_gluon(epsilon) and is_gluon(zeta) and v4 == "f" or v4 == "d"):
    return calculate_6j_with_three_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  #correlates to vertex 1 being triple gluon
  #if(is_gluon(beta) and is_gluon(gamma) and is_gluon(delta) and v1 == "f" or v1 == "d"):
  #  return calculate_6j_with_three_gluon(alpha, conjugate(zeta, n), epsilon, delta, gamma, beta, v4, v2, v3, v1, n, debug)
  #correlates to vertex 2 being triple gluon
  #if(is_gluon(alpha) and is_gluon(gamma) and is_gluon(epsilon) and v2 == "f" or v2 == "d"):
  #  return calculate_6j_with_three_gluon(zeta, beta, delta, gamma, conjugate(epsilon, n), alpha, v1, v4, v3, v2, n, debug)
  #correlates to vertex 3 being triple gluon
  #if(is_gluon(alpha) and is_gluon(beta) and is_gluon(zeta) and v3 == "f" or v3 == "d"):
  #  return calculate_6j_with_three_gluon(conjugate(epsilon, n), conjugate(delta, n), gamma, conjugate(beta, n), conjugate(alpha, n), conjugate(zeta, n), v1, v2, v4, v3, n, debug)
  return False

def calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False, include_coefficient=True):
  """
  Calculates a case 1 6j-symbol.
  """
  if(debug):
    print(debug * "  "+"QuarkGluon", alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4)
  lambdaks = find_intermediate_diagrams(alpha,gamma, False)
  if(len(lambdaks) == 0 or type(lambdaks) == bool):
    return 0 
  if type(v2) == str:
    possible_intermediates = find_intermediate_diagrams(alpha,gamma,False)
    index = int(v2[0])-1
    intermediate = possible_intermediates[index]
    alpha_prime = conjugate(intermediate,n)
    intermediates = find_child_diagram(alpha,gamma)
    alpha_tilde = intermediates[index]
    if alpha_prime in intermediates:
      alpha_tilde = alpha_prime
    if debug:
      print(debug * "  "+"alpha_tilde:",alpha_tilde)
    length = len(intermediates)-1
    if length == 0:
      length = 1
    sign = check_vertex_sign([1,4],alpha,beta,[1,0,0],gamma,alpha_tilde,[1,0,0],1,1,1,1)
    if include_coefficient:
      return calculate_coefficient(alpha,gamma,index,length,n)*Fraction(1,n**2-1)*(sign*calculate_wigner_6j_two_quarks(alpha,beta,[1,0,0],gamma,alpha_tilde,[1,0,0],n=n,debug=debug) - Fraction(kronecker_delta(alpha,gamma),n*get_dimension(alpha,n)))
    if debug:
      print(debug * "  "+"Calculation of QuarkGluon",include_coefficient,Fraction(1,n**2-1),(sign*calculate_wigner_6j_two_quarks(alpha,beta,[1,0,0],gamma,alpha_tilde,[1,0,0],n=n,debug=False) , Fraction(kronecker_delta(alpha,gamma),n*get_dimension(alpha,n))))
    return Fraction(1,n**2-1)*(sign*calculate_wigner_6j_two_quarks(alpha,beta,[1,0,0],gamma,alpha_tilde,[1,0,0],n=n,debug=debug) - Fraction(kronecker_delta(alpha,gamma),n*get_dimension(alpha,n)))
  else:
    intermediates = find_intermediate_diagrams(alpha,gamma,False)
    length = len(intermediates)-1
    if length == 0:
      length = 1
    index = v2-1
    if debug:
      print(debug * "  "+"Calculation of QuarkGluon",include_coefficient,calculate_coefficient(alpha,gamma,length,v2,n),"*",Fraction(1,n**2-1),"*",(Fraction(kronecker_delta(beta,intermediates[index]),get_dimension(beta,n)) ,"-", Fraction(kronecker_delta(alpha,gamma),n*get_dimension(alpha,n))))
    if include_coefficient:
      return calculate_coefficient(alpha,gamma,length,v2,n)*Fraction(1,n**2-1)*(Fraction(kronecker_delta(beta,intermediates[index]),get_dimension(beta,n)) - Fraction(kronecker_delta(alpha,gamma),n*get_dimension(alpha,n)))
    return Fraction(1,n**2-1)*(Fraction(kronecker_delta(beta,intermediates[index]),get_dimension(beta,n)) - Fraction(kronecker_delta(alpha,gamma),n*get_dimension(alpha,n)))

def calculate_6j_with_quark_gluon_conjugate_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Calculates a case 1 6j-symbol.
  """
  if(debug):
    print(debug * "  "+"QuarkGluon", alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4)
  lambdaks = find_intermediate_diagrams(alpha,gamma, False)
  if(len(lambdaks) == 0 or type(lambdaks) == bool):
    return 0 
  if type(v2) == str:
    intermediate = find_child_diagram(alpha,gamma)[int(v2[0])-1]
    return calculate_coefficient(alpha,gamma,int(v2[0]),int(v2[0]),n)*Fraction(1,n**2-1)*(Fraction(kronecker_delta(beta,intermediate),get_dimension(beta,n)) - Fraction(kronecker_delta(alpha,gamma),n*get_dimension(alpha,n)))
  else:
    intermediate = find_intermediate_diagrams(alpha,gamma,False)[v2-1]
    return calculate_coefficient(alpha,gamma,v2,v2,n)*Fraction(1,n**2-1)*(calculate_wigner_6j_two_quarks(alpha,intermediate,[1,0,0],gamma,beta,[1,0,0],n=n) - Fraction(kronecker_delta(alpha,gamma),n*get_dimension(alpha,n)))

def calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False, include_coefficient=True):
  """
  Calculates a case 2 6j-symbol.
  """
  if(debug):
    print((debug-1) * "  "+"QuarkGluonOpposing", alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4) 
  alpha_primes = find_intermediate_diagrams(delta, epsilon)
  
  #print(alpha_primes, delta, epsilon,v4)
  if type(v4) == str:
    intermediates = find_child_diagram(delta,epsilon)
    if debug:
      print(debug * "  "+"Possible alpha_tildes",intermediates)
    index = int(v4[0])
    length = len(intermediates)-1
    if length == 0:
      length = 1
    alpha_prime = intermediates[index-1]
    intermediates_1 = young_multiplication_quark(alpha_prime)
    intermediates_2 = young_multiplication_antiquark(alpha)
    intermediates_3 = young_multiplication_antiquark(beta)
    solution = Fraction()
    if debug:
      print(debug * "  "+"All intermediates, rho should be in all of these",intermediates_1, intermediates_2, intermediates_3)
    for rho in intermediates_1:
      if rho in intermediates_2 and rho in intermediates_3:
        if debug:
          print(debug * "  "+"Rho: ",rho)
        product = Fraction(1,get_dimension(rho,n))
        if debug:
          print(debug * "  "+"Pre-Factor of Triple Product",product)
        sign = check_vertex_sign([1,2,3,4],epsilon,alpha,[1,0,0],rho,alpha_prime,[1,0,0],1,1,1,1)
        product *= sign*calculate_wigner_6j_two_quarks(epsilon,alpha,[1,0,0],rho,alpha_prime,[1,0,0],n=n,debug=debug+1)
        if debug:
          print(debug * "  "+"1st Factor and resulting intermediate product:",sign*calculate_wigner_6j_two_quarks(epsilon,alpha,[1,0,0],rho,alpha_prime,[1,0,0],n=n,debug=False),product)
        product *= calculate_wigner_6j_two_quarks(rho,beta,[1,0,0],delta,alpha_prime,[1,0,0],n=n,debug=debug+1)
        if debug:
          print(debug * "  "+"2nd Factor and resulting intermediate product:",calculate_wigner_6j_two_quarks(rho,beta,[1,0,0],delta,alpha_prime,[1,0,0],n=n,debug=False),product)
        sign2 = check_vertex_sign(2,conjugate(beta,n),conjugate(rho,n),conjugate(alpha,n),[1,0,0],[2,1,0],[1,0,0],1,v3,1,1)
        product *= sign2*Wigner6j(conjugate(beta,n),conjugate(rho,n),conjugate(alpha,n),[1,0,0],[2,1,0],[1,0,0],1,v3,1,1,n,debug+1)
        if debug:
          print(debug * "  "+"3rd Factor and resulting intermediate product:","sign=",sign2,sign2*Wigner6j(conjugate(beta,n),conjugate(rho,n),conjugate(alpha,n),[1,0,0],[2,1,0],[1,0,0],1,v3,1,1,n,False),product)
        solution += product
        if debug:
          print(debug * "  "+"Current total:",solution)
  else:
    if debug:
      print(debug * "  "+"Possible alpha_primes",alpha_primes)
    #if len(alpha_primes)<=v4:
    #  return 0
    alpha_prime = alpha_primes[v4-1]
    intermediates_1 = find_intermediate_diagrams(alpha_prime,alpha_prime)
    intermediates_2 = find_intermediate_diagrams(alpha,alpha)
    intermediates_3 = find_intermediate_diagrams(beta,beta)
    print(debug * "  "+"All intermediates, rho should be in all of these",intermediates_1, intermediates_2, intermediates_3)
    solution = Fraction()
    for rho in intermediates_1:
      if rho in intermediates_2 and rho in intermediates_3:
        if debug:
          print(debug * "  "+"Rho: ",rho)
        product = Fraction(1,get_dimension(rho,n))
        if debug:
          print(debug * "  "+"Pre-Factor of Triple Product",product)
        product *= calculate_wigner_6j_two_quarks(alpha,rho,[1,0,0],alpha_prime,epsilon,[1,0,0],n=n,debug=debug+1)
        if debug:
          print(debug * "  "+"1st Factor and resulting intermediate product:",calculate_wigner_6j_two_quarks(alpha,rho,[1,0,0],alpha_prime,epsilon,[1,0,0],n=n,debug=False),product)
        sign = check_vertex_sign([1,2,3,4],rho,beta,[1,0,0],delta,alpha_prime,[1,0,0],1,1,1,1)
        product *= sign * calculate_wigner_6j_two_quarks(alpha_prime,rho,[1,0,0],beta,delta,[1,0,0],n=n,debug=debug+1)
        if debug:
          print(debug * "  "+"2nd Factor and resulting intermediate product:",sign *calculate_wigner_6j_two_quarks(alpha_prime,rho,[1,0,0],beta,delta,[1,0,0],n=n,debug=False),product)
        sign2 = 1#check_vertex_sign(2,beta,rho,alpha,[1,0,0],[2,1,0],[1,0,0])
        product *= Wigner6j(beta,rho,alpha,[1,0,0],[2,1,0],[1,0,0],1,v3,1,1,n=n,debug=debug+1)
        if debug:
          print(debug * "  "+"3rd Factor and resulting intermediate product and current total:",Wigner6j(beta,rho,alpha,[1,0,0],[2,1,0],[1,0,0],1,v3,1,1,n,False),product)
        if solution == 0:
          solution = product
        else:
          solution += product
        if debug:
          print(debug * "  "+"Current total:",solution)
  return solution

def calculate_6j_with_quark_gluon_opposing_conjugate(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Calculates a case 2 6j-symbol.
  """
  if(debug):
    print((debug-1) * "  "+"QuarkGluonOpposingConjugate", alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4) 
  alpha_primes = find_child_diagram(delta, epsilon)
  if debug:
    print(debug * "  "+"Possible alpha_primes",alpha_primes)
  if len(alpha_primes)<v4:
    return 0
  alpha_prime = alpha_primes[v4-1]
  intermediates_1 = find_intermediate_diagrams(alpha_prime,alpha_prime)
  intermediates_2 = young_multiplication_antiquark(alpha)
  intermediates_3 = young_multiplication_antiquark(beta)
  if debug:
    print(debug * "  "+"All intermediates, rho should be in all of these",intermediates_1, intermediates_2, intermediates_3)
  solution = Fraction()
  for rho in intermediates_1:
    if rho in intermediates_2 and rho in intermediates_2:
      if debug:
        print(debug * "  "+"rho:",rho)
      product = Fraction(1,get_dimension(rho,n))
      if debug:
        print(debug * "  "+"Pre-Factor of Triple Product",product)
      sign = check_vertex_sign([1,2,3,4],epsilon,alpha,[1,0,0],rho,alpha_prime,[1,0,0],1,1,1,1)
      product *= sign*calculate_wigner_6j_two_quarks(epsilon,alpha,[1,0,0],rho,alpha_prime,[1,0,0],n=n,debug=debug+1)
      if debug:
        print(debug * "  "+"1st Factor and resulting intermediate product:",sign*calculate_wigner_6j_two_quarks(epsilon,alpha,[1,0,0],rho,alpha_prime,[1,0,0],n=n,debug=False),product)
      product *= calculate_wigner_6j_two_quarks(rho,beta,[1,0,0],delta,alpha_prime,[1,0,0],n=n,debug=debug+1)
      if debug:
        print(debug * "  "+"2nd Factor and resulting intermediate product:",calculate_wigner_6j_two_quarks(rho,beta,[1,0,0],delta,alpha_prime,[1,0,0],n=n,debug=False),product)
      product *= calculate_6j_with_quark_gluon_vertex(conjugate(beta,3),conjugate(rho,3),conjugate(alpha,3),[1,0,0],[2,1,0],[1,0,0],1,v3,1,1,n,debug+1)
      if debug:
        print(debug * "  "+"3rd Factor and resulting intermediate product:",calculate_6j_with_quark_gluon_vertex(conjugate(beta,3),conjugate(rho,3),conjugate(alpha,3),[1,0,0],[2,1,0],[1,0,0],1,v3,1,1,n),product)
      solution += product
  return solution

def calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False, include_coefficient=True):
  """
  Calculates a case 3 6j-symbol.
  """
  if(debug):
    print("TwoGluon", alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4) 
  alpha_primes = find_intermediate_diagrams(delta, epsilon)
  #print(alpha_primes, delta, epsilon,v4)
  #
  # Conjuagte Case
  #
  if type(v4) == str:
    index = int(v4[0])
    intermediates = find_child_diagram(alpha,gamma)
    length = len(intermediates)-1
    if length == 0:
      length = 1
    alpha_prime = intermediates[index-1]
    intermediates_1 = young_multiplication_gluon(alpha_prime, True)
    intermediates_2 = young_multiplication_antiquark(alpha)
    intermediates_3 = young_multiplication_antiquark(beta)
    solution = Fraction()
    for rho in intermediates_1:
      if rho in intermediates_2 and rho in intermediates_3:
        product = Fraction(1,get_dimension(rho,n))
        if rho == alpha_prime and is_real(rho):
          sum_over_vertices = [1,"1+"]
        if rho == alpha_prime and not is_real(rho):
          sum_over_vertices = [1,"1c",2,"2c"]
        if rho != alpha_prime:
          sum_over_vertices = [1]
        factor_1 = 0
        factor_2 = 0
        for vertex in sum_over_vertices:
          factor_1 += calculate_6j_with_quark_gluon_opposing(epsilon,alpha,[1,0,0],rho,alpha_prime,[1,0,0],1,1,v2,vertex,n)
          factor_2 += calculate_6j_with_quark_gluon_opposing(beta,delta,[1,0,0],alpha_prime,rho,[1,0,0],1,1,v1,vertex,n)
        product *= factor_1
        product *= factor_2
        product *= calculate_6j_with_quark_gluon_vertex(conjugate(alpha,n),conjugate(rho,n),conjugate(beta,n),[1,0,0],[2,1,0],[1,0,0],1,v3,1,1,n=n)
        solution += product
  #
  # NORMAL Case
  #
  else:
    if len(alpha_primes)<=v4:
      return 0
    alpha_prime = alpha_primes[v4-1]
    intermediates_1 = young_multiplication_gluon(alpha_prime)
    intermediates_2 = find_intermediate_diagrams(alpha,alpha)
    intermediates_3 = find_intermediate_diagrams(beta,beta)
    solution = Fraction()
    for rho in intermediates_1:
      if rho in intermediates_2 and rho in intermediates_3:
        product = Fraction(1,get_dimension(rho,n))
        if rho == alpha_prime and is_real(rho):
          sum_over_vertices = [1,"1+"]
        if rho == alpha_prime and not is_real(rho):
          sum_over_vertices = [1,"1c",2,"2c"]
        if rho != alpha_prime:
          sum_over_vertices = [1]
        factor_1 = 0
        factor_2 = 0
        for vertex in sum_over_vertices:
          factor_1 += calculate_6j_with_quark_gluon_opposing(epsilon,alpha,[1,0,0],rho,alpha_prime,[2,1,0],1,1,v2,vertex,n=n)
          factor_2 += calculate_6j_with_quark_gluon_opposing(conjugate(delta,n),conjugate(beta,n),[1,0,0],conjugate(rho,n),conjugate(alpha_prime,n),[2,1,0],1,1,v1,vertex,n=n)
        product *= factor_1
        product *= factor_2
        product *= calculate_6j_with_quark_gluon_vertex(beta,rho,alpha,[1,0,0],[2,1,0],[1,0,0],1,v3,1,1,n=n)
        solution += product
  return solution
  if(debug):
    print("TwoGluon", alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4)
  solution = Fraction()
  for j in range(1,v4+1):
    for k in range(1,v3+1):
      c1 = calculate_coefficient(epsilon,delta,v4,j,n, debug)
      c2 = calculate_coefficient(beta, alpha,v3,k,n, debug)
      if(c1 == "INFINITY" or c2 == "INFINITY"):
        return "COMPLEX INFINITY"
      myks = find_intermediate_diagrams(alpha, beta)
      if((k>= len(myks) and len(myks)>1) or len(myks)==0 or k > len(myks)):
        continue
      lambdajs = find_intermediate_diagrams(delta, epsilon)
      if((j>= len(lambdajs) and len(lambdajs)>1) or len(lambdajs)==0 or j > len(lambdajs)):
        continue
      myk = myks[k-1]
      lambdaj = lambdajs[j-1]
      if(myk == False or lambdaj == False):
        continue
      intermediates = find_intermediate_diagrams(lambdaj, myk)
      #print("Intermediates", intermediates)
      wignerSum = 0
      if(len(intermediates) == 3):
        for i in range(1, len(intermediates)):
          w1 = calculate_6j_with_quark_gluon_opposing(lambdaj,myk,[1,0,0],beta,delta,[2,1,0],1,1,i,v1,n,debug)
          #w2 = calculate_6j_with_quark_gluon_opposing(conjugate(alpha,n),conjugate(epsilon,n),[1,0,0],conjugate(lambdaj,n),conjugate(myk,n),[2,1,0],1,1,v2,v1,n,debug)
          w2 = calculate_6j_with_quark_gluon_opposing(myk,lambdaj,[1,0,0],epsilon,alpha,[2,1,0],1,1,i,v2,n,debug)
          if(w1 == "COMPLEX INFINITY" or w2 == "COMPLEX INFINITY"):
            return "COMPLEX INFINITY"
          wignerSum += w1*w2
      else:
        w1 = calculate_6j_with_quark_gluon_opposing(lambdaj,myk,[1,0,0],beta,delta,[2,1,0],1,1,1,v1,n,debug)
        #w2 = calculate_6j_with_quark_gluon_opposing(conjugate(alpha,n),conjugate(epsilon,n),[1,0,0],conjugate(lambdaj,n),conjugate(myk,n),[2,1,0],1,1,v2,v1,n,debug)
        w2 = calculate_6j_with_quark_gluon_opposing(myk,lambdaj,[1,0,0],epsilon,alpha,[2,1,0],1,1,1,v2,n,debug)
        if(w1 == "COMPLEX INFINITY" or w2 == "COMPLEX INFINITY"):
          return "COMPLEX INFINITY"
        wignerSum += w1*w2
      if(c1 == None or c2 == None or w1 == None or w2 == None):
        continue
      bracket = wignerSum
      if(kronecker_delta(alpha,beta)*kronecker_delta(delta, epsilon)==1):
        bracket -= Fraction(1,(n*get_dimension(delta,n)*get_dimension(alpha,n)))
      partialsum = (c1*c2)/(n**2-1) * bracket
      if(debug):
        print("TwoGluon Calculation: ", c1,"*",c2,"/8*",bracket, "partialsum:", partialsum)
      if(solution == None or partialsum == None):
        continue
      solution += partialsum
  if(debug):
    print("TwoGluon: ", solution)
  if(type(solution)!= Fraction and solution != None):
    print("Inaccuracy in ", alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, solution)
  return solution

def calculate_6j_with_three_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Calculates a case 4 6j-symbol.
  """
  if(debug):
    print("ThreeGluon")
  sum = 0
  if(v4=="f"):
    for j in range(1,v3+1):
      for k in range(1,v1+1):
        for l in range(1,v2+1):
          diagrams = []
          diagrams1 = find_child_diagram(beta, alpha)
          diagrams2 = find_child_diagram(gamma, beta)
          diagrams3 = find_child_diagram(alpha, gamma)
          for diagram in diagrams1:
            if(not diagram in diagrams and diagram in diagrams2 and diagram in diagrams3):
              diagrams.append(diagram)
          partialsum = 0
          if debug:
            print("diagrams, alpha, beta, gamma",diagrams, alpha, beta, gamma)
          if(diagrams == False or diagrams == None):
            continue
          myks = find_intermediate_diagrams(beta, gamma)
          if((k>= len(myks) and len(myks)>1) or len(myks)==0 or k > len(myks)):
            continue
          lambdajs = find_intermediate_diagrams(alpha, beta)
          if((j>= len(lambdajs) and len(lambdajs)>1) or len(lambdajs)==0 or j > len(lambdajs)):
            continue
          nyls = find_intermediate_diagrams(gamma, alpha)
          if((l>= len(nyls) and len(nyls)>1) or len(nyls)==0 or l > len(nyls)):
            continue
          myk = myks[k-1]
          lambdaj = lambdajs[j-1]
          nyl = nyls[l-1]
          if(debug):
            print(lambdaj, myk, nyl)
          if(lambdaj == False or myk == False or nyl == False):
            continue
          for diagram in diagrams:
            w1 = calculate_wigner_6j_two_quarks(alpha, lambdaj, [1,0,0], beta, diagram, [1,0,0],1,1,1,1,n,debug)
            w2 = calculate_wigner_6j_two_quarks(beta, myk, [1,0,0], gamma, diagram, [1,0,0],1,1,1,1,n,debug)
            w3 = calculate_wigner_6j_two_quarks(gamma, nyl,[1,0,0], alpha, diagram, [1,0,0],1,1,1,1,n,debug)
            partialsum += get_dimension(diagram, n) * w1 * w2 * w3
            if(debug):
              print("w1: ",w1,"w2: ",w2,"w3: ",w3,"dimension sigma: ",get_dimension(diagram, n),"partialsum: ", partialsum)
          c1 = calculate_coefficient(beta, alpha, v3, j, n, debug)
          c2 = calculate_coefficient(gamma, beta, v1, k, n, debug)
          c3 = calculate_coefficient(alpha, gamma, v2, l, n, debug)
          if(debug):
            print("c",v3,j,": ",c1,"c",v1,k,": ",c2,"c",v2,l,": ",c3, "alpha: ", alpha, "beta: ", beta, "gamma: ", gamma)
          if(c1 == "INFINITY" or c2 == "INFINITY" or c3 == "INFINITY"):
            return "COMPLEX INFINITY"
          if(c1 == None or c2 == None or c3 == None):
            continue
          sum += c1 * c2 * c3 * (partialsum - Fraction(kronecker_delta(lambdaj,myk)*kronecker_delta(lambdaj,nyl),get_dimension(lambdaj,n)**2))
    sum *= Fraction(1,(n*n-1)**2,1,2*n)
  elif(v4=="d"):
    for j in range(1,v3+1):
      for k in range(1,v1+1):
        for l in range(1,v2+1):
          diagrams = []
          diagrams1 = find_child_diagram(beta, alpha)
          diagrams2 = find_child_diagram(gamma, beta)
          diagrams3 = find_child_diagram(alpha, gamma)
          for diagram in diagrams1:
            if(not diagram in diagrams and diagram in diagrams2 and diagram in diagrams3):
              diagrams.append(diagram)
          partialsum = 0
          myks = find_intermediate_diagrams(beta, gamma)
          if((k>= len(myks) and len(myks)>1) or len(myks)==0 or k > len(myks)):
            continue
          lambdajs = find_intermediate_diagrams(alpha, beta)
          if((j>= len(lambdajs) and len(lambdajs)>1) or len(lambdajs)==0 or j > len(lambdajs)):
            continue
          nyls = find_intermediate_diagrams(alpha, gamma)
          if((l>= len(nyls) and len(nyls)>1) or len(nyls)==0 or l > len(nyls)):
            continue
          myk = myks[k-1]
          lambdaj = lambdajs[j-1]
          nyl = nyls[l-1]
          if(lambdaj == False or myk == False or nyl == False):
            continue
          for diagram in diagrams:
            w1 = calculate_wigner_6j_two_quarks(alpha, lambdaj, [1,0,0], beta, diagram, [1,0,0],1,1,1,1,n,debug)
            w2 = calculate_wigner_6j_two_quarks(beta, myk, [1,0,0], gamma, diagram, [1,0,0],1,1,1,1,n,debug)
            w3 = calculate_wigner_6j_two_quarks(gamma, nyl,[1,0,0], alpha, diagram, [1,0,0],1,1,1,1,n,debug)
            partialsum += get_dimension(diagram, n) * w1 * w2 * w3
            if(debug):
              print("w1: ",w1,"w2: ",w2,"w3: ",w3,"dimension sigma: ",get_dimension(diagram, n),"partialsum: ", partialsum)
          deltaJK = kronecker_delta(lambdaj,myk)
          deltaJL = kronecker_delta(lambdaj,nyl)
          deltaKL = kronecker_delta(myk,nyl)
          bracket = partialsum 
          bracket += Fraction(deltaJK*deltaJL,get_dimension(lambdaj,n)**2)
          bracket += Fraction(4*kronecker_delta(alpha, beta)*kronecker_delta(alpha, gamma),n**2 * get_dimension(alpha, n)**2)
          bracket -= Fraction(2,n) * (Fraction(kronecker_delta(alpha, gamma)*deltaJK, get_dimension(alpha, n)*get_dimension(lambdaj, n)) + Fraction(kronecker_delta(alpha, beta)*deltaKL,get_dimension(alpha, n)* get_dimension(myk,n)) + Fraction(kronecker_delta(beta, gamma)*deltaJL,get_dimension(beta, n)*get_dimension(lambdaj,n)))
          c1 = calculate_coefficient(beta, alpha, v3, j, n, debug)
          c2 = calculate_coefficient(gamma, beta, v1, k, n, debug)
          c3 = calculate_coefficient(alpha, gamma, v2, l, n, debug)
          if(debug):
            print("c",v3,j,": ",c1,"c",v1,k,": ",c2,"c",v2,l,": ",c3, "alpha: ", alpha, "beta: ", beta, "gamma: ", gamma)
          if(c1 == "INFINITY" or c2 == "INFINITY" or c3 == "INFINITY"):
            return "COMPLEX INFINITY"
          if(c1 == None or c2 == None or c3 == None):
            continue
          sum += c1 * c2 * c3 * bracket
    sum *= Fraction(1,(n*n-1)**2,n,2*(n**2-4))
  elif(v4==1):
    w1 = Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "f", n)
    w2 = Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "d", n)
    if(w1 == "COMPLEX INFINITY" or w2 == "COMPLEX INFINITY"):
      return "COMPLEX INFINITY"
    return Fraction(-1,1,n+2,2*n) * w1 + Fraction(1,1,n-2,2*n) * w2
  elif(v4==2):
    w1 = Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "f", n)
    w2 = Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "d", n)
    if(w1 == "COMPLEX INFINITY" or w2 == "COMPLEX INFINITY"):
      return "COMPLEX INFINITY"
    return Fraction(-1,1,n-2,2*n) * w1 - Fraction(1,1,n+2,2*n) * w2
  return sum

def Wigner6j(alpha=None, beta=None, gamma=None, delta=None, epsilon=None, zeta=None, v1=None, v2=None, v3=None, v4=None, n=3, debug=False):
  """
  Overall function to compute 6j-symbols. Assumes directions on the arrows of the 6j-symbol. 
  Returns a single value of debug is False. Prints many intermediate steps if debug in enabled.
  The return is of type Fraction as is defined at the start of the file. 
  """
  if(debug):
    print((debug-1) * "  "+"INPUT: " ,alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n)
  if(isinstance(alpha, bool) or isinstance(beta, bool) or isinstance(gamma, bool) or isinstance(delta, bool) or isinstance(epsilon, bool) or isinstance(zeta, bool)):
    return 0
  solution = False
  if(v1 == "f"):
    return Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, 1, v2, v3, v4, n, debug+1) + Fraction(-1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, 2, v2, v3, v4, n, debug+1)
  if(v1 == "d"):
    return Fraction(1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, 1, v2, v3, v4, n, debug+1) + Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, 2, v2, v3, v4, n, debug+1)
  if(v2 == "f"):
    return Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, debug+1) + Fraction(-1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, 2, v3, v4, n, debug+1)
  if(v2 == "d"):
    return Fraction(1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, debug+1) + Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v2, 2, v3, v4, n, debug+1)
  if(v3 == "f"):
    return Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, 1, v4, n, debug+1) + Fraction(-1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, 2, v4, n, debug+1)
  if(v3 == "d"):
    return Fraction(1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, 1, v4, n, debug+1) + Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v2, v2, 2, v4, n, debug+1)
  if(is_6j_with_two_quarks(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug+1)):
    solution = calculate_wigner_6j_two_quarks(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug+1)
  if(solution != False or (solution == 0 and (type(solution) == int or type(solution) == Fraction))):
    return solution
  elif(is_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug+1) ):
    if alpha != gamma:
      solution = calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug+1)
    if alpha == gamma and not is_real(alpha):
      match v2:
        case 1:
          if get_number_of_boxes(alpha) > get_number_of_boxes(conjugate(alpha,n)):
            solution = calculate_coefficient(conjugate(alpha,n),conjugate(gamma,n),1,1,n) * calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, "1c", v3, v4, n, debug+1, include_coefficient=False)
          else:
            solution = calculate_coefficient(alpha,gamma,1,1,n) * calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug+1, include_coefficient=False)
        case 2:
          if get_number_of_boxes(alpha) > get_number_of_boxes(conjugate(alpha,n)):
            print("Test",calculate_coefficient(conjugate(alpha,n),conjugate(gamma,n),2,1,n) , calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, "1c", v3, v4, n, False, include_coefficient=False),calculate_coefficient(conjugate(alpha,n),conjugate(gamma,n),2,2,n) , calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, "2c", v3, v4, n, debug, include_coefficient=False))
            solution = calculate_coefficient(conjugate(alpha,n),conjugate(gamma,n),2,1,n) * calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, "1c", v3, v4, n, debug+1, include_coefficient=False)
            solution += calculate_coefficient(conjugate(alpha,n),conjugate(gamma,n),2,2,n) * calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, "2c", v3, v4, n, debug+1, include_coefficient=False)
          else:
            print(calculate_coefficient(alpha,gamma,2,1,n) * calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, False, include_coefficient=False),calculate_coefficient(alpha,gamma,2,2,n) * calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, 2, v3, v4, n, debug, include_coefficient=False))
            solution =  calculate_coefficient(alpha,gamma,2,1,n) * calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, debug+1, include_coefficient=False)
            solution += calculate_coefficient(alpha,gamma,2,2,n) * calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, 2, v3, v4, n, debug+1, include_coefficient=False)

    if alpha == gamma and is_real(alpha):
      match v2:
        case "1+":
          print(debug * "  "+"Values: ",calculate_normalization_plus(alpha, find_intermediate_diagrams(alpha,alpha)[0], n, False) , (calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, False, include_coefficient=False) , calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, "1f", v3, v4, n, False, include_coefficient=False)))
          solution = calculate_normalization_plus(alpha, find_intermediate_diagrams(alpha,alpha)[0], n, debug+1) * (calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, debug+1, include_coefficient=False) + calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, "1f", v3, v4, n, debug+1, include_coefficient=False))
        case 1:
          solution = calculate_normalization_minus(alpha, find_intermediate_diagrams(alpha,alpha)[0], n, debug+1) * (calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, debug+1, include_coefficient=False) - calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, "1f", v3, v4, n, debug+1, include_coefficient=False))
  if(solution != False or (solution == 0 and (type(solution) == int or type(solution) == Fraction))):
    return solution
  solution = is_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug+1)
  if(solution):
    if epsilon != delta:
      print(debug * "  "+"Coeff",calculate_coefficient(epsilon,delta,1,1,n))
      solution = calculate_coefficient(epsilon,delta,1,1,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug+1)
    if epsilon == delta and not is_real(epsilon):
      match v4:
        case 1:
          if get_number_of_boxes(epsilon) > get_number_of_boxes(conjugate(epsilon,n)):
            solution = calculate_coefficient(conjugate(epsilon,n),conjugate(delta,n),1,1,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "1c", n, debug, include_coefficient=False)
          else:
            solution = calculate_coefficient(epsilon,delta,1,1,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, 1, n, debug, include_coefficient=False)
        case "1c":
          solution = calculate_coefficient(epsilon,delta,1,1,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "1c", n, debug, include_coefficient=False)
        case 2:
          if get_number_of_boxes(epsilon) > get_number_of_boxes(conjugate(epsilon,n)):
            solution = calculate_coefficient(conjugate(epsilon,n),conjugate(delta,n),2,1,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "1c", n, debug, include_coefficient=False)
            solution += calculate_coefficient(conjugate(epsilon,n),conjugate(delta,n),2,2,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "2c", n, debug, include_coefficient=False)
          else:
            solution = calculate_coefficient(epsilon,delta,2,1,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, 1, n, debug, include_coefficient=False)
            solution += calculate_coefficient(epsilon,delta,2,2,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, 2, n, debug, include_coefficient=False)
        case "2c":
          solution = calculate_coefficient(epsilon,delta,2,1,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "1c", n, debug, include_coefficient=False)
          solution += calculate_coefficient(epsilon,delta,2,2,n) * calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "2c", n, debug, include_coefficient=False)
    if epsilon == delta and is_real(epsilon):
      match v4:
        case "1+":
          print(debug * "  "+"Values",calculate_normalization_plus(epsilon, find_intermediate_diagrams(epsilon,epsilon)[0], n, debug),calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, 1, n, debug, include_coefficient=False),calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "1f", n, debug, include_coefficient=False))
          solution = calculate_normalization_plus(epsilon, find_intermediate_diagrams(epsilon,epsilon)[0], n, debug) * (calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, 1, n, debug, include_coefficient=False) + calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "1f", n, debug, include_coefficient=False))
        case 1:
          solution = calculate_normalization_minus(epsilon, find_intermediate_diagrams(epsilon,epsilon)[0], n, debug) * (calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, 1, n, debug, include_coefficient=False) - calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, "1f", n, debug, include_coefficient=False))
    return solution
  solution = is_6j_with_three_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  if(solution != False or (solution == 0 and (type(solution) == int or type(solution) == Fraction))):
    return solution
  #
  #
  elif(is_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug) ):
    if epsilon != delta:
      solution = calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
    if epsilon == delta and not is_real(delta):
      match v2:
        case 1:
          solution = calculate_coefficient(epsilon,delta,1,1,n) * calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug, include_coefficient=False)
        case "1c":
          solution = calculate_coefficient(epsilon,delta,1,1,n) * calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, "1c", v3, v4, n, debug, include_coefficient=False)
        case 2:
          #if(len(find_intermediate_diagrams(alpha,alpha))<3):
           # return 0
          solution =  calculate_coefficient(epsilon,delta,1,1,n) * calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, debug, include_coefficient=False)
          solution += calculate_coefficient(epsilon,delta,1,1,n) * calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, 2, v3, v4, n, debug, include_coefficient=False)
        case "2c":
          solution = calculate_coefficient(epsilon,delta,1,1,n) * calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, "1c", v3, v4, n, debug, include_coefficient=False)
          solution += calculate_coefficient(epsilon,delta,1,1,n) * calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, "2c", v3, v4, n, debug, include_coefficient=False)
    if epsilon == gamma and is_real(epsilon):
      match v2:
        case "1+":
          solution = calculate_normalization_plus(epsilon, find_intermediate_diagrams(epsilon,epsilon)[0], n, debug) * (calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug, include_coefficient=False) + calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, "1f", v3, v4, n, debug, include_coefficient=False))
        case 1:
          solution = calculate_normalization_minus(epsilon, find_intermediate_diagrams(epsilon,epsilon)[0], n, debug) * (calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug, include_coefficient=False) - calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, "1f", v3, v4, n, debug, include_coefficient=False))
  if(solution != False or (solution == 0 and (type(solution) == int or type(solution) == Fraction))):
    return solution

  #
  #
  #

  return solution



