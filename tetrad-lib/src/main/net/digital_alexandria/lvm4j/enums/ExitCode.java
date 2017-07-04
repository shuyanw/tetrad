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

package net.digital_alexandria.lvm4j.enums;

/**
 * Enum to store exit codes
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public enum ExitCode
{
    // EXIT_ERROR is returned when something went wrong
    EXIT_ERROR(-1),
    // EXIT_SUCCESS is returned when all went well
    EXIT_SUCCESS(0);

    // the exit code
    private final int _CODE;

    ExitCode(final int code) { this._CODE = code; }

    /**
     * Getter for the integer value of the exit code enum.
     *
     * @return returns the exit code as integer
     */
    public int code()
    {
        return _CODE;
    }
}
