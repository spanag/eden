%option noyywrap nounput noinput batch debug
%option reentrant bison-bridge bison-locations
%option yylineno

%{

 #include "LEMS_Expr.tab.h"

#include "stdlib.h"

// https://stackoverflow.com/a/22125500
#define YY_USER_ACTION \
    yylloc->first_line = yylloc->last_line; \
    yylloc->first_column = yylloc->last_column; \
    for(int i = 0; yytext[i] != '\0'; i++) { \
        if(yytext[i] == '\n') { \
            yylloc->last_line++; \
            yylloc->last_column = 0; \
        } \
        else { \
            yylloc->last_column++; \
        } \
    }
//  yyerror (yyscanner, "Unknown character");
%}


%%
 /* see also https://github.com/LEMS/expr-parser/blob/master/src/main/antlr4/imports/CommonLexerRules.g4 */ 
[ \n\t\r\f]+  {    /* skip whitespace */    }

 /* Literals */
[0-9]+((\.[0-9]*)?([eE][\-\+]?[0-9]+)?)?    { /* printf("float\n"); */ yylval->fVal = atof(yytext); return VALUE; }


 /* TODO namespace access operators .:[] */

 /* Operations for expressions */
"*"   { return '*'; }
"/"   { return '/'; }
"-"   { return '-'; }
\+   { return '+'; }
"^"   { return '^'; }
"("   { return '('; }
")"   { return ')'; }
"?"   { return '?'; }
":"   { return ':'; }

".leq."  { return LEQ; }
"<="     { return LEQ; }
".geq."  { return GEQ; }
">="     { return GEQ; }
".lt."   { return LT; }
"<"      { return LT; }
".gt."   { return GT; }
">"      { return GT; }
".eq."   { return EQ; }
"=="     { return EQ; }
".neq."  { return NEQ; }
"!="     { return NEQ; }
".and."  { return AND; }
"and"    { return AND; }
"&&"     { return AND; }
".or."   { return OR; }
"or"     { return OR; }
"||"     { return OR; }
".not."  { return NOT; }
"not"    { return NOT; }
"~"      { return NOT; }

"abs"    { return ABS;    }
"sqrt"   { return SQRT;   }
"sin"    { return SIN;    }
"cos"    { return COS;    }
"tan"    { return TAN;    }
"sinh"   { return SINH;   }
"cosh"   { return COSH;   }
"tanh"   { return TANH;   }
"exp"    { return EXP;    }
"log10"  { return LOG10;  }
"log"    { return LN;     }
"ln"     { return LN;     }
"ceil"   { return CEIL;   }
"floor"  { return FLOOR;  }
"int"    { return INT;    }

"random" { return RANDOM; }

 /* TODO handle symbols too ! */
"H"      { return HFUNC;  }

 /* Names */
[a-zA-Z\_][a-zA-Z0-9\_]*        { yylval->symName = strdup(yytext); /* printf("%s\n",yytext); */ return SYMBOL; }

 /* <<EOF>>               return END_OF_FILE; */

 /* Anything else is garbage */
. { printf("bad \"%s\"\n", yytext); return INVALID; }

%%
