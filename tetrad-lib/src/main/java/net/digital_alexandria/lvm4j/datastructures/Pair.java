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
 * Class that holds a pair of two values that both extend Comparable..
 * <p>
 * Pair implements Comparable. If two Pairs are compared, then the FIRST
 * values are compared.
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 *
 * @param <T> some generic extending Comparable
 * @param <U> some generic extending Comparable
 */
public class Pair<T extends Comparable<T>, U extends Comparable<U>>
    implements Comparable<Pair<T, U>>
{
    // first value of pair
    private final T _T;
    // second value of pair
    private final U _U;

    /**
     * Initialize a object of class Pair.
     * <p>
     * Pair stores two generic Comparables.
     *
     * @param t the first value of the pair
     * @param u the second value of the pair
     */
    public Pair(final T t, final U u)
    {
        this._T = t;
        this._U = u;
    }

    /**
     * Getter for the first value.
     *
     * @return returns the first value
     */
    public final T getFirst()
    {
        return _T;
    }

    /**
     * Getter for the second value.
     *
     * @return returns the second value
     */
    public final U getSecond()
    {
        return _U;
    }

    @Override
    public int compareTo(final Pair<T, U> o)
    {
        return this._T.compareTo(o._T);
    }
}
