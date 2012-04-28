#!/usr/bin/env python
#-*- coding:utf-8 -*-
import re
import sys

DEBUG = True
UNKNOWN_ERROR = 1
Na = 6.02214e23

class ParseError(Exception): pass

class PeriodicTable:
	def __init__(self, protonu=1.007276466812, neutronu=1.00866491600, electronu=0.00054857991):
		self.tableBySym = {}
		self.tableByName = {}
		self.tableByIdx = []
		self.Pu, self.Nu, self.Eu = protonu, neutronu, electronu
	#enddef

	def addElement(self, idx, symbol, name, weight=None, charge=None, unstable=None):
		"""
		Add an element. For optional sections, None works like N/A.
		Charge may be a list of possible charges.
		"""
		elem = (idx, symbol, name, weight, charge, unstable)
		self.tableBySym[symbol] = elem
		self.tableByName[name.lower()] = elem
		if idx >= len(self.tableByIdx):
			self.tableByIdx += [None] * (idx-len(self.tableByIdx)+1)
		#endif
		self.tableByIdx[idx] = elem
	#enddef

	def addElements(self, elements):
		for x in elements:
			if type(x) is list: self.addElement(*x)
			elif type(x) is dict: self.addElement(**x)
		#endfor
	#enddef

	def lookup(self, lu):
		if type(lu) is str:
			llu = lu.lower()
			if llu in self.tableByName: return self.tableByName[llu]
			else: return self.tableBySym[lu]
		else: return self.tableByIdx[lu]
	#endif

	def number(self, lu): return self.lookup(lu)[0]
	def symbol(self, lu): return self.lookup(lu)[1]
	def name(self, lu): return self.lookup(lu)[2]
	def avgu(self, lu): return self.lookup(lu)[3]
	def charge(self, lu): return self.lookup(lu)[4]
	def unstable(self, lu): return self.lookup(lu)[5]

	def calcWeightu(self, sym, amu, charge=None):
		protons = self.number(sym)
		neutrons = amu - protons
		electrons = protons - (self.charge(sym) if charge is None else charge)
		return protons * self.Pu + neutrons * self.Nu + electrons * self.Eu
	#enddef
#endclass

class Molecule:
	elem = re.compile(r'\s*(?:\^?%(num)s)?([A-Z][a-z]*)(?:_?%(num)s)?\s*' % {
		"num": r'([1-9][0-9]*)'
	})
	def __init__(self, pt, expr):
		self._pt = pt
		bits = []
		parents = [bits]
		contents = {}
		tokens = Molecule.elem.split(expr)
		doAppend = True
		for i in range(0,len(tokens),4):
			token = tokens[i]
			try:
				mu, sym, amt = tokens[i+1], tokens[i+2], tokens[i+3]
				if amt is None: amt = 1
				else: amt = int(amt)
				if mu is not None: mu = int(mu)
			except IndexError: doAppend = False
			if token == "(":
				newbit = []
				parents[-1].append([newbit, None, None])
				parents.append(newbit)
			elif token and token[0] == ")":
				parents.pop()
				if len(token) > 1:
					try: amt = int(token[1:])
					except ValueError: raise ParseError(") may only be followed by a number.")
				else: amt = 1
				parents[-1][-1][1] = amt
			elif token != "": raise ParseError("Unknown token (%s)" % token)
			if doAppend: parents[-1].append((sym, amt, mu))
		#endfor
		def recurse(lst, mul):
			for x in lst:
				if type(x[0]) is list: recurse(x[0], mul * x[1])
				elif x[0] in contents: contents[x[0]] += x[1] * mul
				else: contents[x[0]] = x[1] * mul
			#endfor
		#enddef
		recurse(bits, 1)
		self.bits = bits
	#enddef

	def composition(self):
		return self.contents
	#enddef

	def weightu(self, search=None, lst=None):
		if lst is None: lst = self.bits
		total, rtotal = 0, 0
		if type(search) is str: search = [search]
		for x in lst:
			if type(x[0]) is list:
				weight, rweight = self.weightu(search, x[0])
				total += weight * x[1]
				rtotal += rweight * x[1]
			elif search and x[0] not in search: continue
			elif x[2] is None:
				rweight = self._pt.avgu(x[0])
				total += rweight * x[1]
				rtotal += rweight * x[1]
			else:
				total += x[2] * x[1]
				rtotal += self._pt.calcWeightu(x[0], x[2]) * x[1]
			#endif
		#endfor
		return int(total), rtotal
	#enddef

	def __iter__(self): return iter(self.bits)
#endclass

class Mixture:
	def __init__(self, pt, expr):
		self._pt = pt
		self.bits = map((lambda(x): Molecule(pt, x)), expr.split("."))
	#enddef
#endclass

font = {
"lines": 5, "base": 3,
"A0":" _____ ","B0":" _____ ","C0":" _____ ","D0":" ____  ","E0":" _____ ",
"A1":"|  _  |","B1":"| __  |","C1":"|     |","D1":"|    \ ","E1":"|   __|",
"A2":"|     |","B2":"| __ -|","C2":"|   --|","D2":"|  |  |","E2":"|   __|",
"A3":"|__|__|","B3":"|_____|","C3":"|_____|","D3":"|____/ ","E3":"|_____|",
"A4":"       ","B4":"       ","C4":"       ","D4":"       ","E4":"       ",
"F0":" _____ ","G0":" _____ ","H0":" _____ ","I0":" _____ ","J0":"    __ ",
"F1":"|   __|","G1":"|   __|","H1":"|  |  |","I1":"|     |","J1":" __|  |",
"F2":"|   __|","G2":"|  |  |","H2":"|     |","I2":"|-   -|","J2":"|  |  |",
"F3":"|__|   ","G3":"|_____|","H3":"|__|__|","I3":"|_____|","J3":"|_____|",
"F4":"       ","G4":"       ","H4":"       ","I4":"       ","J4":"       ",
"K0":" _____ ","L0":" __    ","M0":" _____ ","N0":" _____ ","O0":" _____ ",
"K1":"|  |  |","L1":"|  |   ","M1":"|     |","N1":"|   | |","O1":"|     |",
"K2":"|    -|","L2":"|  |__ ","M2":"| | | |","N2":"| | | |","O2":"|  |  |",
"K3":"|__|__|","L3":"|_____|","M3":"|_|_|_|","N3":"|_|___|","O3":"|_____|",
"K4":"       ","L4":"       ","M4":"       ","N4":"       ","O4":"       ",
"P0":" _____ ","Q0":" _____ ","R0":" _____ ","S0":" _____ ","T0":" _____ ",
"P1":"|  _  |","Q1":"|     |","R1":"| __  |","S1":"|   __|","T1":"|_   _|",
"P2":"|   __|","Q2":"|  |  |","R2":"|    -|","S2":"|__   |","T2":"  | |  ",
"P3":"|__|   ","Q3":"|__  _|","R3":"|__|__|","S3":"|_____|","T3":"  |_|  ",
"P4":"       ","Q4":"   |__|","R4":"       ","S4":"       ","T4":"       ",
"U0":" _____ ","V0":" _____ ","W0":" _ _ _ ","X0":" __ __ ","Y0":" __ __ ",
"U1":"|  |  |","V1":"|  |  |","W1":"| | | |","X1":"|  |  |","Y1":"|  |  |",
"U2":"|  |  |","V2":"|  |  |","W2":"| | | |","X2":"|-   -|","Y2":"|_   _|",
"U3":"|_____|","V3":" \___/ ","W3":"|_____|","X3":"|__|__|","Y3":"  |_|  ",
"U4":"       ","V4":"       ","W4":"       ","X4":"       ","Y4":"       ",
"Z0":" _____ ",
"Z1":"|__   |",
"Z2":"|   __|",
"Z3":"|_____|",
"Z4":"       ",

"a0":"     ","b0":" _   ","c0":"     ","d0":"   _ ","e0":"     ","f0":" ___ ",
"a1":" ___ ","b1":"| |_ ","c1":" ___ ","d1":" _| |","e1":" ___ ","f1":"|  _|",
"a2":"| .'|","b2":"| . |","c2":"|  _|","d2":"| . |","e2":"| -_|","f2":"|  _|",
"a3":"|__,|","b3":"|___|","c3":"|___|","d3":"|___|","e3":"|___|","f3":"|_|  ",
"a4":"     ","b4":"     ","c4":"     ","d4":"     ","e4":"     ","f4":"     ",
"g0":"     ","h0":" _   ","i0":" _ ","j0":"   _ ","k0":" _   ","l0":" _ ",
"g1":" ___ ","h1":"| |_ ","i1":"|_|","j1":"  |_|","k1":"| |_ ","l1":"| |",
"g2":"| . |","h2":"|   |","i2":"| |","j2":"  | |","k2":"| '_|","l2":"| |",
"g3":"|_  |","h3":"|_|_|","i3":"|_|","j3":" _| |","k3":"|_,_|","l3":"|_|",
"g4":"|___|","h4":"     ","i4":"   ","j4":"|___|","k4":"     ","l4":"   ",
"m0":"       ","n0":"     ","o0":"     ","p0":"     ","q0":"     ","r0":"     ",
"m1":" _____ ","n1":" ___ ","o1":" ___ ","p1":" ___ ","q1":" ___ ","r1":" ___ ",
"m2":"|     |","n2":"|   |","o2":"| . |","p2":"| . |","q2":"| . |","r2":"|  _|",
"m3":"|_|_|_|","n3":"|_|_|","o3":"|___|","p3":"|  _|","q3":"|_  |","r3":"|_|  ",
"m4":"       ","n4":"     ","o4":"     ","p4":"|_|  ","q4":"  |_|","r4":"     ",
"s0":"     ","t0":" _   ","u0":"     ","v0":"     ","w0":"       ","x0":"     ",
"s1":" ___ ","t1":"| |_ ","u1":" _ _ ","v1":" _ _ ","w1":" _ _ _ ","x1":" _ _ ",
"s2":"|_ -|","t2":"|  _|","u2":"| | |","v2":"| | |","w2":"| | | |","x2":"|_'_|",
"s3":"|___|","t3":"|_|  ","u3":"|___|","v3":" \_/ ","w3":"|_____|","x3":"|_,_|",
"s4":"     ","t4":"     ","u4":"     ","v4":"     ","w4":"       ","x4":"     ",
"y0":"     ","z0":"     ",
"y1":" _ _ ","z1":" ___ ",
"y2":"| | |","z2":"|- _|",
"y3":"|_  |","z3":"|___|",
"y4":"|___|","z4":"     ",

"_0":"       ","(0":"    /| ",")0":" |\    ",
"_1":"       ","(1":"   / | ",")1":" | \   ",
"_2":"       ","(2":"  / /  ",")2":"  \ \  ",
"_3":"  ___  ","(3":"  \ \  ",")3":"  / /  ",
"_4":" |___| ","(4":"   \_| ",")4":" |_/   ",
}
def spoolIsotopic(sym, num=None, amu=None, amt=None, charge=None, spool=None):
	global font
	if spool is None: spool = [""] * font["lines"]

	# Spool left-side numbers (AMU top, Atomic Number bottom)
	pad = 0
	if num: num = str(num)
	else: num = ""
	if amu:
		amu = str(amu)
		spool[0] += (" " * max(len(num) - len(amu),0)) + amu
		pad = len(amu)
	else:
		amu = ""
		spool[0] += " " * len(amu)
	#endif
	spool[font["base"]] += (" " * max(len(amu) - len(num),0)) + num
	pad = max(pad,len(num))
	if pad:
		pad = " " * pad
		for i in range(1,font["lines"]):
			if i != font["base"]: spool[i] += pad
		#endfor
	#endif

	# Spool symbol
	spoolWord(sym, spool)

	# Spool charge and amount
	pad = 0
	if amt: amt = "" if amt == 1 else str(amt)
	else: amt = ""
	if charge:
		charge = chargeToStr(int(charge))
		spool[0] += charge + (" " * max(len(amt) - len(charge),0))
		pad = len(charge)
	else:
		charge = ""
		spool[0] += " " * len(amt)
	#endif
	spool[font["base"]] += amt + (" " * max(len(charge) - len(amt),0))
	pad = max(pad,len(amt))
	if pad:
		pad = " " * pad
		for i in range(1,font["lines"]):
			if i != font["base"]: spool[i] += pad
		#endfor
	#endif

	return spool
#enddef

def spoolMolecule(compound, spool=None):
	global font
	if spool is None: spool = [""] * font["lines"]
	for x in compound:
		if type(x[0]) is list:
			spoolWord("(", spool)
			spoolMolecule(x[0], spool)
			spoolWord(")", spool)
			if x[1] > 1:
				amt = str(x[1])
				pad = " " * len(amt)
				for i in range(font["lines"]):
					if i != font["base"]: spool[i] += pad
					else: spool[i] += amt
				#endfor
			#endif
		else: spoolIsotopic(x[0], amt=x[1], amu=x[2], spool=spool)
	#endfor
	return spool
#enddef

def spoolWord(word, spool=None):
	global font
	if spool is None: spool = [""] * font["lines"]
	first = True
	for x in word:
		# Note: Underscore must always be defined!!
		try: font[x+"0"]
		except IndexError: x = "_"
		if first: first, collided = False, True
		else:
			# First verify no colisions...
			collided = False
			for i in range(font["lines"]):
				nc = font[x+str(i)][0]
				lc = spool[i][-1]
				if nc != lc and nc != " " and lc != " ":
					collided = True
					break
				#endif
			#endfor
		#endif

		if collided:
			# Now spool
			for i in range(font["lines"]): spool[i] += font[x+str(i)]
		else:
			# Or smush and spool
			for i in range(font["lines"]):
				line = font[x+str(i)]
				if spool[i][-1] == " ": spool[i] = spool[i][0:-1] + line[0]
				spool[i] += line[1:]
			#endfor
		#endif
	#endfor
	return spool
#enddef

def spoolSpace(space=5, spool=None):
	global font
	if spool is None: spool = [""] * font["lines"]
	space = " " * space
	map((lambda(x): x + space), spool)
	return spool
#enddef

def printSpool(spool): print("\n".join(spool))

def parseCharge(charge):
	if   charge == "":      charge =  0
	elif charge == "-":     charge = -1
	elif charge == "+":     charge =  1
	elif charge[-1] == "-": charge = -int(charge[0:-1])
	elif charge[-1] == "+": charge =  int(charge[0:-1])
	else: charge = int(charge)
	return charge
#enddef

def chargeToStr(charge):
	if   charge == 0:  charge = ""
	elif charge == 1:  charge = "+"
	elif charge == -1: charge = "-"
	elif charge > 1:   charge = "%i+" % charge
	elif charge < -1:  charge = "%i-" % abs(charge)
	return charge
#enddef

################################################################################
elementTable = PeriodicTable()
elementTable.addElements([
[1,"H","Hydrogen",1.00794,1,False],[2,"He","Helium",4.002602,0,False],
[3,"Li","Lithium",6.941,1,False],[4,"Be","Beryllium",9.012182,2,False],
[5,"B","Boron",10.811,0,False],[6,"C","Carbon",12.0107,0,False],
[7,"N","Nitrogen",14.0067,-3,False],[8,"O","Oxygen",15.9994,-2,False],
[9,"F","Fluorine",18.9984032,-1,False],[10,"Ne","Neon",20.1797,0,False],
[11,"Na","Sodium",22.98976928,1,False],[12,"Mg","Magnesium",24.305,2,False],
[13,"Al","Aluminium",26.9815386,3,False],[14,"Si","Silicon",28.0855,0,False],
[15,"P","Phosphorus",30.973762,-3,False],[16,"S","Sulfur",32.065,-2,False],
[17,"Cl","Chlorine",35.453,-1,False],[18,"Ar","Argon",39.948,0,False],
[19,"K","Potassium",39.0983,1,False],[20,"Ca","Calcium",40.078,2,False],
[21,"Sc","Scandium",44.955912,3,False],[22,"Ti","Titanium",47.867,[4,3],False],
[23,"V","Vanadium",50.9415,[3,5],False],[24,"Cr","Chromium",51.9961,[3,2],False],
[25,"Mn","Manganese",54.938045,[2,4],False],[26,"Fe","Iron",55.845,[3,2],False],
[27,"Co","Cobalt",58.933195,[2,3],False],[28,"Ni","Nickel",58.6934,[2,3],False],
[29,"Cu","Copper",63.546,[2,1],False],[30,"Zn","Zinc",65.38,2,False],
[31,"Ga","Gallium",69.723,3,False],[32,"Ge","Germanium",72.63,4,False],
[33,"As","Arsenic",74.9216,-3,False],[34,"Se","Selenium",78.96,-2,False],
[35,"Br","Bromine",79.904,-1,False],[36,"Kr","Krypton",83.798,0,False],
[37,"Rb","Rubidium",85.4678,1,False],[38,"Sr","Strontium",87.62,2,False],
[39,"Y","Yttrium",88.90585,3,False],[40,"Zr","Zirconium",91.224,4,False],
[41,"Nb","Niobium",92.90638,[5,3],False],[42,"Mo","Molybdenum",95.96,6,False],
[43,"Tc","Technetium",98,7,True],[44,"Ru","Ruthenium",101.07,[3,4],False],
[45,"Rh","Rhodium",102.9055,3,False],[46,"Pd","Palladium",106.42,[2,4],False],
[47,"Ag","Silver",107.8682,1,False],[48,"Cd","Cadmium",112.411,2,False],
[49,"In","Indium",114.818,3,False],[50,"Sn","Tin",118.71,[4,2],False],
[51,"Sb","Antimony",121.76,[3,5],False],[52,"Te","Tellurium",127.6,-2,False],
[53,"I","Iodine",126.90447,-1,False],[54,"Xe","Xenon",131.293,0,False],
[55,"Cs","Caesium",132.9054519,1,False],[56,"Ba","Barium",137.327,2,False],
[57,"La","Lanthanum",138.90547,3,False],[58,"Ce","Cerium",140.116,3,False],
[59,"Pr","Praseodymium",140.90765,3,False],[60,"Nd","Neodymium",144.242,3,False],
[61,"Pm","Promethium",145,3,True],[62,"Sm","Samarium",150.36,[3,2],False],
[63,"Eu","Europium",151.964,[3,2],False],[64,"Gd","Gadolinium",157.25,3,False],
[65,"Tb","Terbium",158.92535,3,False],[66,"Dy","Dysprosium",162.5,3,False],
[67,"Ho","Holmium",164.93032,3,False],[68,"Er","Erbium",167.259,3,False],
[69,"Tm","Thulium",168.93421,3,False],[70,"Yb","Ytterbium",173.054,[3,2],False],
[71,"Lu","Lutetium",174.9668,3,False],[72,"Hf","Hafnium",178.49,4,False],
[73,"Ta","Tantalum",180.94788,5,False],[74,"W","Tungsten",183.84,6,False],
[75,"Re","Rhenium",186.207,7,False],[76,"Os","Osmium",190.23,4,False],
[77,"Ir","Iridium",192.217,4,False],[78,"Pt","Platinum",195.084,[4,2],False],
[79,"Au","Gold",196.966569,[3,1],False],[80,"Hg","Mercury",200.59,[2,1],False],
[81,"Tl","Thallium",204.3833,[1,3],False],[82,"Pb","Lead",207.2,[2,4],False],
[83,"Bi","Bismuth",208.9804,[3,5],False],[84,"Po","Polonium",209,[2,4],True],
[85,"At","Astatine",210,-1,True],[86,"Rn","Radon",222,0,True],
[87,"Fr","Francium",223,1,True],[88,"Ra","Radium",226,2,True],
[89,"Ac","Actinium",227,3,True],[90,"Th","Thorium",232.03806,4,False],
[91,"Pa","Protactinium",231.03588,[5,4],False],[92,"U","Uranium",238.02891,[6,4],False],
[93,"Np","Neptunium",237,5,True],[94,"Pu","Plutonium",244,[4,6],True],
[95,"Am","Americium",243,[3,4],True],[96,"Cm","Curium",247,3,True],
[97,"Bk","Berkelium",247,[3,4],True],[98,"Cf","Californium",251,3,True],
[99,"Es","Einsteinium",252,3,True],[100,"Fm","Fermium",257,3,True],
[101,"Md","Mendelevium",258,[2,3],True],[102,"No","Nobelium",259,[2,3],True],
[103,"Lr","Lawrencium",262,3,True],[104,"Rf","Rutherfordium",267,None,True],
[105,"Db","Dubnium",268,None,True],[106,"Sg","Seaborgium",271,None,True],
[107,"Bh","Bohrium",272,None,True],[108,"Hs","Hassium",270,None,True],
[109,"Mt","Meitnerium",276,None,True],[110,"Ds","Darmstadtium",281,None,True],
[111,"Rg","Roentgenium",280,None,True],[112,"Cn","Copernicium",285,None,True],
[113,"Uut","Ununtrium",284,None,True],[114,"Fl","Flerovium",289,None,True],
[115,"Uup","Ununpentium",288,None,True],[116,"Lv","Livermorium",293,None,True],
[117,"Uus","Ununseptium",294,None,True],[118,"Uuo","Ununoctium",294,None,True],
])

def cast_input(prompt, blank=None, cast=int, noblank=False):
	while 1:
		tmp = raw_input(prompt)
		if tmp == "" and not noblank: return blank
		try: return cast(tmp)
		except ValueError: pass
	#endwhile
#enddef

################################################################################
def main(args, argc):
	global elementTable
	table = elementTable
	post = False
	while 1:
		if post: print("")
		else: post = True
		print("====== Main Menu ======"
		"\n1) Print isotopic notation"
		"\n2) Enter compound"
		"\n3) Lookup"
		"\nq) Quit"
		)
		choice = raw_input("Choice: ")
		if choice == "1": # Print isotopic notation
			symbol = raw_input("Symbol [skip]: ")
			try:
				number = table.number(symbol)
				tmpnum = cast_input("Number (/ for unshown) [%s]: " % number,
					cast=(lambda(x): x if x == "/" else int(x))
				)
				if tmpnum == "/": number = None
				elif tmpnum is not None: number = tmpnum
			except KeyError:
				number = cast_input("Number [%s]: " % (
					"cancel" if symbol == "" else "unshown"
				))
			#endtry
			if symbol == "":
				if number: symbol = table.symbol(number)
				else: continue
			#endif
			amu = cast_input("Weight (AMU) [unshown]: ")
			# Ensure proper output format
			charge = parseCharge(raw_input("Charge [0]: "))
			amount = cast_input("Amount [1]: ")
			if amount <= 1: amount = None
			print("")
			printSpool(spoolIsotopic(symbol, number, amu, amount, charge))
		elif choice == "2": # Enter compound
			while 1:
				toparse = raw_input("Enter compound (? for help): ")
				if toparse == "?":
					print("For each element you may enter (in order):"
					"\n ^# (Optional) - Weight (AMU)"
					"\n Sy - Element Symbol"
					"\n # (Optional) - Amount represented"
					"\nYou may also group elements with parenthesis. A number "
					"may follow the closing parenthesis. Do not use spaces."
					"\nType q to quit."
					)
				elif toparse == "q": break
				elif toparse != "":
					compound = Molecule(table, toparse)
					break
				#endif
			#endwhile
			if toparse == "q": continue
			while 1:
				print("\n===== Compound Menu ====="
				"\n1) Display"
				"\n2) Weight"
				"\n3) I have X atoms..."
				"\n4) I have X moles..."
				"\n5) I have X grams..."
				"\nq) Quit"
				)
				choice = raw_input("Choice: ")
				if choice == "1": printSpool(spoolMolecule(compound))
				elif choice == "2":
					of = raw_input("Of? [all]: ")
					if of == "": of = None
					else: of = re.findall(r'[A-Z][a-z]*', of)
					w = compound.weightu(of)
					print("Simple: %s AMU\nPrecise: %s AMU" % w)
					if of: print("Mass percent: %s%%" % (w[1] / compound.weightu()[1] * 100))
				elif choice == "3":
					atoms = cast_input("Atoms (e notation): ", cast=float)
					if atoms is None: continue
					print("You have %e moles." % (atoms / Na))
					print("You have %e grams." % (atoms * compound.weightu()[1] / Na))
				elif choice == "4":
					moles = cast_input("Moles (e notation): ", cast=float)
					if moles is None: continue
					print("You have %e atoms." % (moles * Na))
					print("You have %e grams." % (moles * compound.weightu()[1]))
				elif choice == "5":
					grams = cast_input("Grams (e notation): ", cast=float)
					if grams is None: continue
					amu = compound.weightu()[1]
					print("You have %e atoms." % (grams * Na / amu))
					print("You have %e moles." % (grams / amu))
				elif choice == "q": break
			#endwhile
		elif choice == "3":
			try:
				idx, symbol, name, weight, charge, unstable = (
				table.lookup(raw_input("Atomic Number, Symbol, or Name: ")))
				if type(charge) is not list: charge = [charge]
				if unstable: weight = "(%s)" % weight
				print("===== %s =====" % name)
				for x in charge:
					printSpool(spoolIsotopic(symbol, idx, weight, None, x))
				#endfor
			except KeyError, IndexError: continue
		elif choice == "q": return 0
	#endwhile
#enddef

if __name__ == "__main__":
	try: sys.exit(main(sys.argv, len(sys.argv)))
	except KeyboardInterrupt:
		print("")
		sys.exit(0)
	except Exception:
		if DEBUG: raise
		else: sys.exit(UNKNOWN_ERROR)
	#endtry
#endif
