/**
 * lvm4j: a Java implementation of various latent variable models.
 *
 * Copyright (C) 2015 - 2016 Simon Dirmeier
 *
 * This file is part of lvm4j.
 * <p>
 * lvm4j is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * lvm4j is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with lvm4j.  If not, see <http://www.gnu.org/licenses/>.
 */


package net.digital_alexandria.lvm4j.datastructures;


/**
 * Class that holds a triple of three values. The first ones are comparable
 * while the third one is not.
 * <p>
 * Triple is a subclass of Pair, so can be compared.
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 *
 * @param <T> some generic extending Comparable
 * @param <U> some generic extending Comparable
 * @param <V> some generic
 */
public final class Triple<T extends Comparable<T>, U extends
    Comparable<U>, V> extends Pair<T, U>
{
    // third value of triple
    private final V _V;

    /**
     * Initialize a object of class Triple.
     * <p>
     * Triple stores two generic Comparables and one random generic.
     *
     * @param t the first value of the triple
     * @param u the second value of the triple
     * @param v the third value of the tripl
     */
    public Triple(final T t, final U u, final V v)
    {
        super(t, u);
        this._V = v;
    }

    /**
     * Getter for the third value.
     *
     * @return returns the third value
     */
    public V getThird()
    {
        return _V;
    }
}
