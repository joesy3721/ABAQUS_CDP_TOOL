# ABAQUS_CDP_TOOL
This is a useful MatLab based tool to easily set up the Concrete Damaged Plasticity Model available in Abaqus. 

Check out the tutorial on Youtube:
https://www.youtube.com/watch?v=NtLSh5OsB5c

Place the .fig and the .m file in the same folder and open and run the .m file in MatLab. 

Start with specifying the name for your material in accordance to your section definition in Abaqus.
Simply enter material parameters like compressive and tensile strength, Young's modulus and tensile Mode I fracture energy and determine your desired damage limit. 
With the help of the Abaqus Manual you should be able to understand all parameters.
Since the compressive behavior is not regularized in this Concrete Model (in contrast to the tensile behavior) you can enter a characteristic element length to controll mesh sensitivity. Therefore you must determine the ratio of compressive to tensile fracture energy which usually lies between 250 and 350.
You can determine the compressive damage - inelastic strain curve by entering your desired damage levels at the compressive strength and at the strain limit according to EC1992-1-1 or Model Code 2012. Keep in mind that damage curves must be monotonously increasing.
You can determine the tensile damage - inelastic displacement curve by entering your desired damage limit.
You can export your data so you can simply include it in your ABAQUS .inp file by:

** MATERIALS
**
*INCLUDE, Input = yourexportfilename.txt

References:
Compressive behavior:
EN 1992-1-1:2004 + AC:2008 + AC:2010 + A1:2014 (D), Eq.(3.14), Tab.3.1
Tensile behavior:
G. Hofstetter, G. Meschke, Numerical Modeling of Concrete Cracking, Springer-Verlag Wien, 2011 (Doi: 10.1007/978-3-7091-0897-0)
Compressive/Tensile Fracture Energy ratio:
Hikaru Nakamura, Takeshi Higai, Compressive Fracture Energy and Fracture Zone Length of Concrete, 2011
Abaqus Manual:
https://help.3ds.com/HelpProductsDS.aspx
