import math
import numpy as np


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
      return (s.nominator == o.nominator and s.denominator == o.denominator and s.nominatorRoot == o.nominatorRoot and s.denominatorRoot == o.denominatorRoot)
    if isinstance(other, (float,int)):
      return self.denominator * other * self.denominatorRoot == self.nominator *math.sqrt(self.nominatorRoot)

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

def get_partitions(numbers, sum, max_number=None):
  """
  Generates all possible partitions with a certain amount of numbers and given total sum. max_number can be provided to limit the maximum possible digit.
  This is used to generate Young diagrams.
  """
  if sum == 0:
    return [np.zeros(numbers, dtype=int)]
  elif numbers == 1:
    return [[sum]]
  start = sum
  if max_number:
    start = min(sum, max_number)
  partitions = False
  for i in range(start, math.ceil(sum / numbers)-1, -1):
    if type(partitions) == bool:
      partitions = np.insert(get_partitions(numbers - 1, sum - i, i), 0, i, axis=1)
    else:
      partitions = np.append(partitions,np.insert(get_partitions(numbers - 1, sum - i, i), 0, i, axis=1),axis=0)
  return partitions

def create_young_tableaus(n, to_sum = 10, max_amount = 1000, includeSinglet = False):
  """
  Creates all Young Tableaus with an amount of boxes up to n + to_sum and return them in order in a list.
  """
  tableaus = False
  sum = 1
  if(includeSinglet):
    sum = 0
  created = 0
  while sum <= to_sum and created <= max_amount:
    if type(tableaus) == bool:
      tableaus = get_partitions(n, sum)
    else:
      tableaus = np.append(tableaus, get_partitions(n, sum), axis=0)
      sum += 1
    created += 1
  return [tableau.tolist() for tableau in tableaus[1:]]

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
  if(alpha == [0,0,0]): #TODO Compatible with n
    return 0
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

def young_addition(alpha, beta):
  """
  Adds the diagram obtained by addind rowlengths of both provided diagrams.
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
    if(diag != False):
      diagrams.append(diag)
  return diagrams

def young_multiplication_gluon(alpha):
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
      if(j != False):
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

def calculate_wigner_6j_two_quarks(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Calculates a case 0 6j-symbol.
  """
  if(debug):
    print("TwoQuarkSymbol", alpha, beta, gamma, delta, epsilon, zeta)
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
    return calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  #if(find_changed_index(beta, alpha) != -1 and is_gluon(delta) and is_quark(epsilon) and is_quark(conjugate(zeta, n)) and v2 == 1 and v3 == 1 and v4 == 1):
  #  return calculate_6j_with_quark_gluon_vertex(gamma, alpha, beta, conjugate(zeta, n), delta, epsilon, v3, v1, v2, v4, n, debug)
  #TODO there are more cases
  if(is_gluon(epsilon) and is_quark(conjugate(delta, n)) and is_quark(conjugate(zeta, n)) and v1 == 1 and v2 == 1 and v4 == 1):
    return calculate_6j_with_quark_gluon_vertex(gamma, beta, alpha, [1,0,0], epsilon, [1,0,0], v3, v2, v1, v4, n, debug)
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
  if(is_quark(gamma) and is_gluon(zeta) and v1 == 1 and v2 == 1 and (v4 == 1 or v4 == 2)):
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

def calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Calculates a case 1 6j-symbol.
  """
  if(debug):
    print("QuarkGluon", alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4)
  lambdaks = find_intermediate_diagrams(alpha,gamma, False)
  if(len(lambdaks) == 0 or type(lambdaks) == bool):
    return 0 
  if(alpha != gamma and beta == lambdaks[0]):
    return Fraction(1,1,1,get_dimension(beta, n)*(n**2-1)).reduce()
  elif(alpha != gamma):
    return 0
  intermediates = find_intermediate_diagrams(alpha, alpha, False)
  k=0
  for i in range(len(intermediates)):
    if intermediates[i] == beta:
      k = i+1
  if(k==0):
    print("Problem in Quark Gluon", alpha, beta, gamma, intermediates)
    return 0
  #k = find_changed_index(alpha, beta)+1
  if(v2==1):
    if debug:
      print(k, scalar_product(1,k,alpha,n),scalar_product(1,1,alpha,n))
    return scalar_product(1,k,alpha,n) / Fraction.sqrt(scalar_product(1,1,alpha,n))
  if(v2==2):
    #print(Fraction.sqrt(scalar_product(1,1,alpha,n)/(scalar_product(1,1,alpha,n)*scalar_product(2,2,alpha,n)-Fraction.sq(scalar_product(1,2,alpha,n)))))
    #print((scalar_product(2,k,alpha,n)),(scalar_product(1,2,alpha,n)/scalar_product(1,1,alpha,n)*scalar_product(1,k,alpha,n)))
    divisor = scalar_product(1,1,alpha,n)*scalar_product(2,2,alpha,n) - Fraction.sq(scalar_product(1,2,alpha,n))
    bracket = (scalar_product(2,k,alpha,n)-((scalar_product(1,2,alpha,n)/scalar_product(1,1,alpha,n))*scalar_product(1,k,alpha,n)))
    if(debug):
      print("s11:",scalar_product(1,1,alpha,n),"s22",scalar_product(2,2,alpha,n) ,"s12",scalar_product(1,2,alpha,n),"s12",scalar_product(2,1,alpha,n),"k", k, "divisor: ", divisor, "bracket",bracket)

    if(bracket == 0):
      return 0
    if(divisor == 0):
      return "COMPLEX INFINITY"
    return Fraction.sqrt(scalar_product(1,1,alpha,n)/divisor)*bracket

def calculate_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Calculates a case 2 6j-symbol.
  """
  if(debug):
    print("QuarkGluonOpposing", alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4) 
  solution = Fraction()
  for j in range(1,v4+1):
    for k in range(1,v3+1):
      c1 = calculate_coefficient(epsilon,delta,v4,j,n, debug)
      c2 = calculate_coefficient(beta, alpha,v3,k,n, debug)
      if(c1 == "INFINITY" or c2 == "INFINITY"):
        return "COMPLEX INFINITY"
      myks = find_intermediate_diagrams(alpha, beta)
      #print(myks)
      if((k>= len(myks) and len(myks)>1) or len(myks)==0 or k > len(myks)):
        continue
      lambdajs = find_intermediate_diagrams(delta, epsilon)
      #print(lambdajs)
      if((j>= len(lambdajs) and len(lambdajs)>1) or len(lambdajs)==0 or j > len(lambdajs)):
        continue
      myk = myks[k-1]
      lambdaj = lambdajs[j-1]
      #if(debug):
        #print("myk", myk,"lambdaj", lambdaj)
      if(myk == False or lambdaj == False):
        continue
      w1 = calculate_wigner_6j_two_quarks(conjugate(lambdaj, n),conjugate(delta, n),[1,0,0],conjugate(beta, n),conjugate(myk, n),[1,0,0],1,1,1,1,n,debug)
      w2 = calculate_wigner_6j_two_quarks(alpha,myk,[1,0,0],lambdaj,epsilon,[1,0,0],1,1,1,1,n,debug)
      #print(w1,"test",w2)
      if(w1 == "COMPLEX INFINITY" or w2 == "COMPLEX INFINITY"):
        return "COMPLEX INFINITY"
      if(c1 == None or c2 == None or w1 == None or w2 == None):
        continue
      bracket = w1*w2
      partialsum = (c1*c2)/(n**2-1) * bracket
      if(kronecker_delta(alpha,beta)*kronecker_delta(delta, epsilon)==1):
        bracket -= Fraction(1,(n*get_dimension(delta,n)*get_dimension(alpha,n)))
      if(solution == None or partialsum == None):
        continue
      solution += partialsum
      if(debug):
        if(kronecker_delta(alpha,beta)*kronecker_delta(delta, epsilon)==1):
          print("QuarkGluon Calculation: ", partialsum, "=", c1,"*",c2,"/",(n**2-1)," *",bracket, ", last term is: ",w1,"*",w2,"-",Fraction(1,(n*get_dimension(delta,n)*get_dimension(alpha,n))))
        else:
          print("QuarkGluon Calculation: ", partialsum, "=", c1,"*",c2,"/",(n**2-1)," *",bracket, ", last term is: ",w1,"*",w2)
  if(debug):
    print("QuarkGluon6j: ",solution)
  return solution

def calculate_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug=False):
  """
  Calculates a case 3 6j-symbol.
  """
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
    print("INPUT: " ,alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n)
  if(isinstance(alpha, bool) or isinstance(beta, bool) or isinstance(gamma, bool) or isinstance(delta, bool) or isinstance(epsilon, bool) or isinstance(zeta, bool)):
    return 0
  solution = False
  if(v1 == "f"):
    return Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, 1, v2, v3, v4, n, debug) + Fraction(-1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, 2, v2, v3, v4, n, debug)
  if(v1 == "d"):
    return Fraction(1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, 1, v2, v3, v4, n, debug) + Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, 2, v2, v3, v4, n, debug)
  if(v2 == "f"):
    return Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, debug) + Fraction(-1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, 2, v3, v4, n, debug)
  if(v2 == "d"):
    return Fraction(1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, 1, v3, v4, n, debug) + Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v2, 2, v3, v4, n, debug)
  if(v3 == "f"):
    return Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, 1, v4, n, debug) + Fraction(-1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, 2, v4, n, debug)
  if(v3 == "d"):
    return Fraction(1,1,n-2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, 1, v4, n, debug) + Fraction(-1,1,n+2,2*n) * Wigner6j(alpha, beta, gamma, delta, epsilon, zeta, v2, v2, 2, v4, n, debug)

  if(is_6j_with_two_quarks(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)):
    solution = calculate_wigner_6j_two_quarks(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  if(solution != False or (solution == 0 and (type(solution) == int or type(solution) == Fraction))):
    return solution
  elif(is_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug) ):
    solution = calculate_6j_with_quark_gluon_vertex(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  if(solution != False or (solution == 0 and (type(solution) == int or type(solution) == Fraction))):
    return solution
  solution = is_6j_with_quark_gluon_opposing(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  if(solution != False or (solution == 0 and (type(solution) == int or type(solution) == Fraction))):
    return solution
  solution = is_6j_with_three_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  if(solution != False or (solution == 0 and (type(solution) == int or type(solution) == Fraction))):
    return solution
  solution = is_6j_with_two_gluon(alpha, beta, gamma, delta, epsilon, zeta, v1, v2, v3, v4, n, debug)
  if(solution != False or (solution == 0 and (type(solution) == int or type(solution) == Fraction))):
    return solution
  return solution





