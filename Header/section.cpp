//
// Created by 조준홍 on 2022/02/18.
//

#include "section.h"

ALG::section::section() = default;

ALG::section::section(std::string name) : name(std::move(name)) {}

const std::string &ALG::section::getName() const {
    return name;
}

void ALG::section::setName(const std::string &string) {
    section::name = string;
}