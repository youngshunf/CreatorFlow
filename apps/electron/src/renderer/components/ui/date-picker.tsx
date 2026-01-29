/**
 * DatePicker - Date selection component with calendar popover.
 *
 * Combines Calendar and Popover for a user-friendly date selection experience.
 */

import * as React from 'react'
import { format, parse, isValid } from 'date-fns'
import { zhCN } from 'date-fns/locale'
import { CalendarIcon } from 'lucide-react'

import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import { Calendar } from '@/components/ui/calendar'
import {
  Popover,
  PopoverContent,
  PopoverTrigger,
} from '@/components/ui/popover'

export interface DatePickerProps {
  /** Selected date value (YYYY-MM-DD format) */
  value?: string
  /** Callback when date changes */
  onChange?: (value: string) => void
  /** Placeholder text when no date selected */
  placeholder?: string
  /** Disabled state */
  disabled?: boolean
  /** Additional className for the trigger button */
  className?: string
  /** Date format for display (default: yyyy年MM月dd日) */
  displayFormat?: string
  /** Minimum selectable date */
  minDate?: Date
  /** Maximum selectable date */
  maxDate?: Date
}

/**
 * DatePicker - Single date selection with calendar popover
 *
 * @example
 * <DatePicker
 *   value={birthday}
 *   onChange={setBirthday}
 *   placeholder="选择日期"
 * />
 */
export function DatePicker({
  value,
  onChange,
  placeholder = '选择日期',
  disabled = false,
  className,
  displayFormat = 'yyyy年MM月dd日',
  minDate,
  maxDate,
}: DatePickerProps) {
  const [open, setOpen] = React.useState(false)

  // Parse string value to Date object
  const dateValue = React.useMemo(() => {
    if (!value) return undefined
    const parsed = parse(value, 'yyyy-MM-dd', new Date())
    return isValid(parsed) ? parsed : undefined
  }, [value])

  // Handle date selection
  const handleSelect = React.useCallback(
    (date: Date | undefined) => {
      if (date && onChange) {
        onChange(format(date, 'yyyy-MM-dd'))
      }
      setOpen(false)
    },
    [onChange]
  )

  return (
    <Popover open={open} onOpenChange={setOpen}>
      <PopoverTrigger asChild>
        <Button
          variant="outline"
          disabled={disabled}
          className={cn(
            'w-full justify-start text-left font-normal',
            !dateValue && 'text-muted-foreground',
            className
          )}
        >
          <CalendarIcon className="mr-2 size-4" />
          {dateValue ? format(dateValue, displayFormat, { locale: zhCN }) : placeholder}
        </Button>
      </PopoverTrigger>
      <PopoverContent className="w-[280px] p-0" align="start">
        <Calendar
          mode="single"
          selected={dateValue}
          onSelect={handleSelect}
          captionLayout="dropdown"
          defaultMonth={dateValue}
          disabled={(date) => {
            if (minDate && date < minDate) return true
            if (maxDate && date > maxDate) return true
            return false
          }}
          fromYear={1900}
          toYear={new Date().getFullYear()}
        />
      </PopoverContent>
    </Popover>
  )
}
