//
// Created by gregor on 25.06.19.
//

#include "Argument.h"

/*
 * specialization for TP=std:string
 *
 */
template<>
void Argument<std::string>::treat_default_sval() {
    std::string default_sval = this->data_[0];
    std::vector<std::string> snippets;

    size_t found = default_sval.find(" ");
    while (found != std::string::npos) {
        snippets.push_back(default_sval.substr(0, found));
        default_sval = default_sval.substr(found + 1);
        found = default_sval.find(" ");
    }
    snippets.push_back(default_sval);

    this->nargs_ = snippets.size();
    delete[] this->data_;
    this->data_ = new std::string[this->nargs_];
    for (int i = 0; i < std::max<int>(this->nargs_, 1); i++) {
        this->data_[i] = snippets[i];
    } /* i */
}

/*
 * \brief specialization for string
 *
 */
template<>
int Argument<std::string>::parse(int argc, char *argv[]) {
    bool custom(false);
    if (this->pos_ == 0) {
        std::stringstream inss(std::stringstream::in | std::stringstream::out);

        for (int i = 0; i < argc; i++) {
            std::string arg;
            inss.clear();
            inss << argv[i];
            inss >> arg;
            if (this->is(arg)) {
                custom = true;
                if (this->nargs_ == -1) { count_nargs_(i, argc, argv); }
                /* assert that option argument exists */
                assert(i + this->nargs_ < argc);
                for (int j = 0; j < this->nargs_; j++) {
                    inss.clear();
                    inss << argv[i + j + 1];
                    inss >> this->data_[j];
                } /* j */
            }
        } /* i */
        if (!custom) { treat_default_sval(); }
    } else if (argc > this->pos_) {
        std::stringstream inss(std::stringstream::in | std::stringstream::out);
        inss << argv[this->pos_];
        inss >> this->data_[0];
    }

    parsed = true;
    return 0;
}


/*
 * \brief specialization for bool
 */
template<>
int Argument<bool>::parse(int argc, char *argv[]) {
    if (this->pos_ == 0) {
        std::stringstream inss(std::stringstream::in | std::stringstream::out);

        for (int i = 0; i < argc; i++) {
            std::string arg;
            inss.clear();
            inss << argv[i];
            inss >> arg;
            if (this->is(arg)) {
                if (this->nargs_ > 0) {
                    assert(i + this->nargs_ < argc);
                    for (int j = 0; j < this->nargs_; j++) {
                        inss.clear();
                        inss << argv[i + j + 1];
                        inss >> this->data_[j];
                    }
                } else {
                    this->data_[0] = true;
                }
            }
        } /* i */
    } else if (argc > this->pos_) {
        std::stringstream inss(std::stringstream::in | std::stringstream::out);
        inss << argv[this->pos_];
        inss >> this->data_[0];
    }
    return 0;
}