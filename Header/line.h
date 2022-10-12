//
// Created by 조준홍 on 2022/02/18.
//

#ifndef ALG2D_LINE_H
#define ALG2D_LINE_H

#include "matrix.h"

namespace ALG {
    class line2D {
    private:
        vector start{}, end{}, normal{}, tangent{};
        double length{};
        bool iscloseline{false};

    public:
        line2D();

        line2D(const vector& start, const vector &anEnd);

        line2D(double s0, double s1, double e0, double e1);

        void calcProperties();

        [[nodiscard]] vector &getStart();

        [[nodiscard]] const vector &getStart() const;

        void setStart(const vector &vector);

        [[nodiscard]] vector &getAnEnd();

        [[nodiscard]] const vector &getAnEnd() const;

        void setAnEnd(const vector &vector);

        [[nodiscard]] vector &getNormal();

        [[nodiscard]] const vector &getNormal() const;

        [[nodiscard]] vector &getTangent();

        [[nodiscard]] const vector &getTangent() const;

        [[nodiscard]] double getLength() const;

        [[nodiscard]] bool getIscloseLine() const;

        bool iscross(const line2D &src, ALG::vector &vec);

        virtual ~line2D();
    };

}


#endif //ALG2D_LINE_H
