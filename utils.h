//
// Created by federico on 18/07/17.
//

#ifndef DISCRETELOGARITHMS_UTILS_H
#define DISCRETELOGARITHMS_UTILS_H

struct element {
    mpz_t value;
    int exp;
    struct element *next;
    struct element *prev;
};

struct list {
    struct element *HEAD;
    struct element *TAIL;
    int count;
};

void addelem(mpz_t value, struct list *list);
void freelist(struct list *list);
void printlist(struct list *list);
void returnlist(struct list *linkedlist, mpz_t *returnedlist);
struct list *init_list();
struct element *addnewfactor(mpz_t value, struct list *list);

#endif //DISCRETELOGARITHMS_UTILS_H

