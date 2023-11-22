#include "LEMS_Expr.h"

#include "LEMS_Expr.tab.h"
#include "LEMS_Expr.yy.h"


void TermTable::printTree(const TermTable &tab, int node, int tabDepth){
	auto &term = tab[node];
	for(int i = 0; i < tabDepth; i++) printf("    ");
	printf("< %d >\n",node);
	for(int i = 0; i < tabDepth; i++) printf("    ");
	switch(term.type){
	case Term::VALUE: {
		printf("%f\n",term.value);
	break;}
	case Term::SYMBOL: {
		printf("#%d ", term.symbol);
		printf("%s\n", tab.symbol_refs[term.symbol].c_str());
	break;}
	case Term::PLUS: {
		printf("+ \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::MINUS: {
		printf("- \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::TIMES: {
		printf("* \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::DIVIDE: {
		printf("/ \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::POWER: {
		printf("^ \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	
	case Term::LEQ: {
		printf("<= \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::GEQ: {
		printf(">= \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::LT: {
		printf("< \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::GT: {
		printf("> \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::EQ: {
		printf("== \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::NEQ: {
		printf("!= \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::AND: {
		printf("and \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::OR: {
		printf("or \n");
		printTree(tab, term.left, tabDepth + 1);
		printTree(tab, term.right, tabDepth + 1);
	break;}
	
	case Term::UMINUS: {
		printf("- \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::UPLUS: {
		printf("+ \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	
	case Term::NOT: {
		printf("not \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	
	case Term::ABS: {
		printf("abs \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::SQRT: {
		printf("sqrt \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::SIN: {
		printf("sin \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::COS: {
		printf("cos \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::TAN: {
		printf("tan \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::SINH: {
		printf("sinh \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::COSH: {
		printf("cosh \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::TANH: {
		printf("tanh \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::EXP: {
		printf("exp \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::LOG10: {
		printf("log10 \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::LN: {
		printf("ln \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::CEIL: {
		printf("ceil \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::FLOOR: {
		printf("floor \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::RANDOM: {
		printf("random \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::HFUNC: {
		printf("H \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	case Term::INT: {
		printf("int \n");
		printTree(tab, term.right, tabDepth + 1);
	break;}
	
	default: {
		printf("unknown term!!!\n");
	break;}
	}
}

bool ParseLemsExpression( const char *sInput, TermTable &tab ){
	yyscan_t scanner;	
	if(yylex_init(&scanner)){
		return false;
	}
	
	YY_BUFFER_STATE buf = yy_scan_string(sInput, scanner);
	
	if(yyparse( scanner, &tab, sInput ) ){
		return false;
	}
	yy_delete_buffer(buf, scanner);
	yylex_destroy(scanner);
	
	return true;
}
