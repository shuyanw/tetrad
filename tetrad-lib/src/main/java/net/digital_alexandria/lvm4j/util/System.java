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


import net.digital_alexandria.lvm4j.enums.ExitCode;

/**
 * Class that is used for exiting, printing to stderr, etc.
 *
 * @author Simon Dirmeier {@literal s@simon-dirmeier.net}
 */
public final class System
{
    /** private constructor to avoid instantiation **/
    private System() {}

    /**
     * Print a message to the error stream and exit the program with a
     * specific error code.
     *
     * @param message the message to be printed
     * @param exitCode the error code that is returned
     */
    public static void exit(final String message, final ExitCode exitCode)
    {
        java.lang.System.err.println(message);
        java.lang.System.exit(exitCode.code());
    }
}
