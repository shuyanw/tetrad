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


package net.digital_alexandria.lvm4j.util;

import java.util.ArrayList;
import java.util.List;


/**
 * Class Combinatorial can be used to create all combinations of letters, etc.
 * for a given length.
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public final class Combinatorial
{
    /** private constructor to avoid instantiation **/
    private Combinatorial() {}

    /**
     * Get every combination of strings of lengths for a given array of chars
     * .<br>
     * E.g: for array [a,b] and length 2 create a List of strings : [aa, ab,
     * ba, bb].
     *
     * @param array the chars that should be combinatorically created.
     * @param length the length of every string that is returned
     *
     * @return returns a list of strings
     */
    public static List<String> combinatorial(final char[] array,
                                             final int length)
    {
        List<String> list = new ArrayList<>();
        combinatorial(new StringBuilder(), list, array, length);
        return list;
    }

    /**
     * Appender for the combinatorial expansion of a char array
     *
     * @param builder current combination as string builder
     * @param list list with all the combinations
     * @param array array of elements that are combined
     * @param length maximal length of single combination
     */
    private static void combinatorial(final StringBuilder builder,
                                      final List<String> list,
                                      final char[] array, final int length)
    {
        if (builder.length() != 0 && builder.length() < length)
            list.add(builder.toString());
        if (builder.length() == length)
            list.add(builder.toString());
        else
        {
            for (char ac : array)
            {
                combinatorial(new StringBuilder(builder).append(ac),
                              list, array, length);
            }
        }
    }
}
