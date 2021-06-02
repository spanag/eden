#ifndef LEMS_EXPR_H
#define LEMS_EXPR_H

#include <string>
#include <vector>
#include <algorithm>

struct Term{

	enum Type{
		NONE,
		VALUE,
		SYMBOL,
		
		//binary operators
		PLUS,
		MINUS,
		TIMES,
		DIVIDE,
		POWER,
		
		LEQ,
		GEQ,
		LT,
		GT,
		EQ,
		NEQ,
		AND,
		OR,
		
		// unary functions
		UMINUS,
		UPLUS,
		
		NOT,
		
		ABS   ,
		SQRT  ,
		SIN   ,
		COS   ,
		TAN   ,
		SINH  ,
		COSH  ,
		TANH  ,
		EXP   ,
		LOG10 ,
		LN    ,
		CEIL  ,
		FLOOR ,
		RANDOM,
		HFUNC ,
			
	}type;
	
		
	int left;
	int right;
	double value;
	int symbol;
	

	
	Term(Type _type, int _l, int _r){
		type = _type;
		left = _l;
		right = _r;
	}
	Term(Type _type, int _r){
		type = _type;
		right = _r;
	}
	Term(double _v){
		type = VALUE;
		value = _v;
	}
	static Term Symbol(int _s){
		Term term;
		term.type = SYMBOL;
		term.symbol = _s;
		return term;
	}
	
	bool isBinaryOperator() const {
		return false
			|| type == PLUS
			|| type == MINUS
			|| type == TIMES 
			|| type == DIVIDE 
			|| type == POWER

			|| type == LEQ
			|| type == GEQ
			|| type == LT
			|| type == GT 
			|| type == EQ 
			|| type == NEQ 
			|| type == AND 
			|| type == OR 
		;
	}
	bool isUnary() const {
		return ( isUnaryOperator() || isUnaryFunction() );
	}
	bool isUnaryOperator() const {
		return false
			|| type == UMINUS
			|| type == UPLUS
			|| type == NOT
		;
	}
	bool isUnaryFunction() const {
		return false
			|| type == ABS   
			|| type == SQRT  
			|| type == SIN   
			|| type == COS   
			|| type == TAN   
			|| type == SINH  
			|| type == COSH  
			|| type == TANH  
			|| type == EXP   
			|| type == LN    
			|| type == LOG10 
			|| type == CEIL  
			|| type == FLOOR 
			|| type == RANDOM
			|| type == HFUNC 
		;
	}
	
	bool isBoolean() const {
		return false
			|| type == LEQ
			|| type == GEQ
			|| type == LT 
			|| type == GT 
			|| type == EQ 
			|| type == NEQ
			|| type == AND
			|| type == OR 
			|| type == NOT 
		;
	}
	bool isNumeric() const {
		return !( isBoolean() || isUndefined() ); // for now
	}
	bool isUndefined() const {
		return type == NONE;
	}
private:
	Term(){} //internal use only
};
struct TermTable{
	std::vector<Term> terms;
	int expression_root;
	
	std::vector<std::string> symbol_refs;
	
	int add( const Term &new_term ){
		int new_id = terms.size();
		//printf("term %d -> %d \n", new_term.type, new_id);
		terms.push_back(new_term);		
		return new_id;
	};
	
	int addSymbol( const char * sSymName ){
		std::string symname = sSymName;
		int symid;
		auto it = std::find(symbol_refs.begin(), symbol_refs.end(), symname);
		if( it != symbol_refs.end() ){
			symid = it - symbol_refs.begin();
		}
		else{
			symid = symbol_refs.size();
			symbol_refs.push_back(symname);
		}
		
		int new_id = terms.size();
		terms.push_back( Term::Symbol(symid) );		
		return new_id;
	}

	Term &operator[](int index){
		return terms.at(index);
	}
	const Term &operator[](int index) const {
		return terms.at(index);
	}
	
	TermTable(){
		expression_root = -1;
	}
	
	//for user diagnostics
	// TODO print infix
	
	// for debug
	static void printTree(const TermTable &tab, int node, int tabDepth = 0);
};

bool ParseLemsExpression( const char *sInput, TermTable &tab );

#endif
