//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_SECTION_H
#define ALG2D_SECTION_H


#include "lineOthers.h"

namespace ALG {
    class section : public std::vector<boundaryLine2D> {
    private:
        std::string name{};
    public:
        section();

        explicit section(std::string name);

        [[nodiscard]] const std::string &getName() const;

        void setName(const std::string &string);
    };

}


#endif //ALG2D_SECTION_H
