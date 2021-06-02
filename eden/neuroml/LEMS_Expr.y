%require "3.0"
%locations

// error reporting luxury
%define parse.trace
%define parse.error verbose

%define api.pure full
// could use %param too, for the args used in both
%param   { yyscan_t scanner }
%lex-param  { TermTable *tab } { const char *original_string } 
%parse-param { TermTable *tab } { const char *original_string } 

%code requires{

#define YYDEBUG 1

// This is the best place to write dependency code required for YYSTYPE and YYLTYPE. In other words, it’s the best place to define types referenced in %union directives. If you use #define to override Bison’s default YYSTYPE and YYLTYPE definitions, then it is also the best place. However you should rather %define api.value.type and api.location.type. 
#include "LEMS_Expr.h"

// TODO
// #define YYSTYPE val_ptr
// #define YYSTYPE val_ptr // term ID

// just in case, define yyscan_t to use it
#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
typedef void *yyscan_t;
#endif
	
}

%code provides{

#define YY_DECL int yylex(YYSTYPE * yylval_param, YYLTYPE * yylloc_param , yyscan_t yyscanner, struct TermTable *tab, const char *original_string)
	YY_DECL;

}

%code {
//#include "flex.h"

// error handler TODO
void yyerror(YYLTYPE *yylloc, yyscan_t scann, TermTable *tab, const char *original_string,  const char *msg) {
	int fc = yylloc->first_column - 1;
	int lc = yylloc->last_column - 1;
	printf("Error on byte %d ~ %d: %s \n", fc + 1, lc + 1 -1, msg);
	// also show string
	printf("%s\n", original_string);
	for(int i = 0; i < fc; i++ ) printf(" ");
	if( fc <= lc){
		for(int i = fc; i < lc ;i++) printf("^");
	}
	printf("\n");
}


}

%union {
	int termId;					 /* semantic term, not a raw token */
    double fVal;                 /* floating-point value */
    const char *symName;         /* symbol allocated with strdup */
};

%token <fVal> VALUE "number"
%token <symName> SYMBOL "symbol"

%token LEQ GEQ LT GT EQ NEQ AND OR NOT
%token ABS SQRT SIN COS TAN SINH COSH TANH EXP LOG10 LN CEIL FLOOR RANDOM HFUNC

%token NUM INVALID

%token END_OF_FILE "end of file"
%token END 0 "end of expression"

%type <termId> anyexpr "arithmetic expression"
%type <termId> numexpr "numeric expression"
%type <termId> num_function "built-in function"
%type <termId> binexpr "boolean expression"

%type <termId> symexpr "symbol expression"

// Operator associations and precedence.

%left OR
%left AND

%left EQ NEQ
%left GEQ LEQ

%left '+' '-'
%left '*' '/'
%right '^'



%nonassoc ABS SQRT SIN COS TAN SINH COSH TANH EXP LOG10 LN CEIL FLOOR RANDOM HFUNC
%right NOT
%right UMINUS "unary -"

%%
 // see also https://github.com/LEMS/expr-parser/blob/master/src/main/antlr4/parser/LEMSExpression.g4
anyexpr
	: numexpr END { tab->expression_root = $1; }
	| binexpr END { tab->expression_root = $1; } 
	;

numexpr
	: numexpr '+' numexpr { /* printf("plus %d %d\n", $1, $3); */ $$ = tab->add(Term( Term::PLUS  , $1, $3 )); }
	| numexpr '-' numexpr { $$ = tab->add(Term( Term::MINUS , $1, $3 )); }
	| numexpr '*' numexpr { $$ = tab->add(Term( Term::TIMES , $1, $3 )); }
	| numexpr '/' numexpr { $$ = tab->add(Term( Term::DIVIDE, $1, $3 )); }

	| numexpr '^' numexpr { $$ = tab->add(Term( Term::POWER, $1, $3 )); }
	// TODO ternary
	// | binexpr '?' numexpr ':' numexpr { $$ = tab->add(Term( Term::TERNARY, $1, $3, $5 )); }
	| '-' numexpr %prec UMINUS { $$ = tab->add(Term( Term::UMINUS, $2 )); }
	| '+' numexpr %prec UMINUS { $$ = tab->add(Term( Term::UPLUS, $2 )); }
	| num_function    { $$ = $1; }
	| '(' numexpr ')' { $$ = $2; }
	| symexpr { $$ = $1; }
	| VALUE  { $$ = tab->add( Term($1) ); }
	; 

binexpr
	: binexpr AND binexpr { $$ = tab->add(Term( Term::AND , $1, $3 )); }
	| binexpr OR  binexpr { $$ = tab->add(Term( Term::OR  , $1, $3 )); }
	|         NOT binexpr { $$ = tab->add(Term( Term::NOT ,     $2 )); }
	| numexpr LEQ numexpr { $$ = tab->add(Term( Term::LEQ , $1, $3 )); }
	| numexpr GEQ numexpr { $$ = tab->add(Term( Term::GEQ , $1, $3 )); }
	| numexpr LT  numexpr { $$ = tab->add(Term( Term::LT  , $1, $3 )); }
	| numexpr GT  numexpr { $$ = tab->add(Term( Term::GT  , $1, $3 )); }
	| numexpr EQ  numexpr { $$ = tab->add(Term( Term::EQ  , $1, $3 )); }
	| numexpr NEQ numexpr { $$ = tab->add(Term( Term::NEQ , $1, $3 )); }
	| '(' binexpr ')' { $$ = $2; }	
	// | symexpr { $$ = $1; } a single symbol would be ambiguous as to whether boolean or numeric!
	;

num_function
	: ABS    '(' numexpr ')' { $$ = tab->add(Term( Term::ABS   , $3 )); }
	| SQRT   '(' numexpr ')' { $$ = tab->add(Term( Term::SQRT  , $3 )); }
	| SIN    '(' numexpr ')' { $$ = tab->add(Term( Term::SIN   , $3 )); }
	| COS    '(' numexpr ')' { $$ = tab->add(Term( Term::COS   , $3 )); }
	| TAN    '(' numexpr ')' { $$ = tab->add(Term( Term::TAN   , $3 )); }
	| SINH   '(' numexpr ')' { $$ = tab->add(Term( Term::SINH  , $3 )); }
	| COSH   '(' numexpr ')' { $$ = tab->add(Term( Term::COSH  , $3 )); }
	| TANH   '(' numexpr ')' { $$ = tab->add(Term( Term::TANH  , $3 )); }
	| EXP    '(' numexpr ')' { $$ = tab->add(Term( Term::EXP   , $3 )); }
	| LOG10  '(' numexpr ')' { $$ = tab->add(Term( Term::LOG10 , $3 )); }
	| LN     '(' numexpr ')' { $$ = tab->add(Term( Term::LN    , $3 )); }
	| CEIL   '(' numexpr ')' { $$ = tab->add(Term( Term::CEIL  , $3 )); }
	| FLOOR  '(' numexpr ')' { $$ = tab->add(Term( Term::FLOOR , $3 )); }
	| RANDOM '(' numexpr ')' { $$ = tab->add(Term( Term::RANDOM, $3 )); }
	| HFUNC  '(' numexpr ')' { $$ = tab->add(Term( Term::HFUNC , $3 )); }
	;

symexpr
	: SYMBOL { $$ = tab->addSymbol( $1 ); free( (void *) $1); } 
	;
%%
