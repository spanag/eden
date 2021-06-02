#include <stdio.h>
#include "LEMS_Expr.h"

int main(int argc, char *argv[]) {
	if(argc < 2){
		printf("need to put an expression string\n");
		return 6;
	}
	
	TermTable tab;
	
	const char *sInput = argv[1];
	
	if( !ParseLemsExpression( sInput, tab ) ){
		return -1;
	}
  
	if(tab.expression_root < 0){
		printf("could not find expression root");
		return -1;
	}
	
	printf("start %d symbols %zd\n", tab.expression_root, tab.symbol_refs.size());
	TermTable::printTree(tab, tab.expression_root, 0);
	/*
	for(int i = 0; i < tab.terms.size(); i++){
		printf("tab %d\n", i);
		TermTable::printTree(tab, i, 0);
	}*/
	
	return 0;
}
